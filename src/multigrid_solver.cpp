#define _USE_MATH_DEFINES
#include "gravomg/multigrid_solver.h"
#include "gravomg/utility.h"

#include <cmath>
#include <numeric>
#include <chrono>

#include <Eigen/Eigenvalues>

namespace GravoMG {
	/* Constructor */
	MultigridSolver::MultigridSolver(
		Eigen::MatrixXd& V, Eigen::MatrixXi& neigh, Eigen::SparseMatrix<double>& M) : 
		V(V), neigh(neigh), M(M) {
		Minv = (Eigen::SparseMatrix<double>) M.diagonal().cwiseInverse().asDiagonal();

		V0 = V;
		hierarchyTiming["n_vertices"] = V.rows();
	}

	/* ====================== BASIC FUNCTIONS  ======================= */
	MultigridSolver::~MultigridSolver()
	{
		// std::cout << "Deleting data!\n";
		V.resize(0, 0);
		V0.resize(0, 0);
		neigh.resize(0, 0);
		DoF.clear(); DoF.shrink_to_fit();
		U.clear(); U.shrink_to_fit();
		Abar.clear(); Abar.shrink_to_fit();
		Mbar.clear(); Mbar.shrink_to_fit();
		neighHierarchy.clear(); neighHierarchy.shrink_to_fit();
		nearestSource.clear(); nearestSource.shrink_to_fit();
		geodDistance.clear(); geodDistance.shrink_to_fit();
		samples.clear(); samples.shrink_to_fit();
		PointsToSampleMap.clear(); PointsToSampleMap.shrink_to_fit();
	}

	// Build the hierarchy
	void MultigridSolver::buildHierarchy()
	{
		nearestSource.clear();
		nearestSource.shrink_to_fit();

		plf::nanotimer timer;
		timer.start();
		constructProlongation();
		hierarchyTiming["hierarchy"] = timer.get_elapsed_ms();
	}

	void MultigridSolver::constructProlongation()
	{
		// Create prolongation operator for level k+1 to k
		// Points in current level
		Eigen::MatrixXd levelPoints = V;
		Eigen::MatrixXd levelNormals = normals;
		// Points in coarser level
		Eigen::MatrixXd samplePoints;
		// Neighborhood data structure
		Eigen::MatrixXi neighLevelK = neigh;

		// Compute initial radius
		double densityRatio = std::sqrt(ratio);
		int nLevel1 = int(levelPoints.rows() / ratio);

		std::random_device					rd;
		std::default_random_engine			generator(rd());

		hierarchyTiming["sampling"] = 0.0;
		hierarchyTiming["cluster"] = 0.0;
		hierarchyTiming["next_neighborhood"] = 0.0;
		hierarchyTiming["next_positions"] = 0.0;
		hierarchyTiming["triangle_finding"] = 0.0;
		hierarchyTiming["triangle_selection"] = 0.0;
		hierarchyTiming["levels"] = DoF.size() - 1.0;
		// For each level
		int k = 0;
		DoF.clear();
		DoF.shrink_to_fit();
		DoF.push_back(levelPoints.rows());
		while (levelPoints.rows() > lowBound && k < 10) {
			double radius = std::cbrt(ratio) * computeAverageEdgeLength(levelPoints, neighLevelK);
			// -- Setting up variables
			// Data structure for neighbors inside level k
			std::vector<std::set<int>> neighborsList;

			// List of triplets to build prolongation operator U
			std::vector<Eigen::Triplet<double>> AllTriplet, UNeighAllTriplet;

			// The nearest coarse point for each fine point and the distances computed
			Eigen::VectorXd D(levelPoints.rows());
			D.setConstant(std::numeric_limits<double>::max());
			nearestSource.push_back(std::vector<size_t>(levelPoints.rows()));

			// -- Sample a subset for level k + 1
			if (verbose) printf("__Constructing Prolongation Operator for level = %d using closest triangle. \n", k);

			// Sample points that will be part of the coarser level with Poisson disk sampling
			if (verbose) std::cout << "Obtaining subset from the finer level\n";
			
			plf::nanotimer samplingTimer;
			samplingTimer.start();

			switch (samplingStrategy) {
				case FASTDISK:
					samples.push_back(fastDiskSample(levelPoints, neighLevelK, radius, D, nearestSource[k]));
					DoF.push_back(samples[k].size());
					break;
				case RANDOM:
					DoF.push_back(DoF[k] / ratio);
					samples.push_back(std::vector<int>(DoF[k]));
					std::iota(samples[k].begin(), samples[k].end(), 0);
					std::shuffle(samples[k].begin(), samples[k].end(), generator);
					samples[k].resize(DoF[k + 1]);
					break;
				case MIS:
					samples.push_back(maximumDeltaIndependentSet(levelPoints, neighLevelK, radius, D, nearestSource[k]));
					DoF.push_back(samples[k].size());
					break;
			}

			if (samples[k].size() < lowBound) {
				nearestSource.pop_back();
				break;
			}

			if (verbose) cout << "Actual number found: " << samples[k].size() << endl;
			DoF[k + 1] = samples[k].size();

			hierarchyTiming["sampling"] += samplingTimer.get_elapsed_ms();

			plf::nanotimer clusterTimer;
			clusterTimer.start();			

			// -- Compute distance from fine points to coarse points and get closest coarse point
			// using distances from MIS if computed before
			constructDijkstraWithCluster(levelPoints, samples[k], neighLevelK, k, D, nearestSource[k]); // Stores result in nearestSource[k]
			hierarchyTiming["cluster"] += clusterTimer.get_elapsed_ms();

			plf::nanotimer nextNeighTimer;
			nextNeighTimer.start();

			// Create neighborhood for the next level
			neighborsList.resize(DoF[k + 1]);
			for (int fineIdx = 0; fineIdx < DoF[k]; ++fineIdx) {
				for (int j = 0; j < neighLevelK.cols(); ++j) {
					int neighIdx = neighLevelK(fineIdx, j);
					if (neighIdx < 0) break;
					if (nearestSource[k][fineIdx] != nearestSource[k][neighIdx]) {
						neighborsList[nearestSource[k][fineIdx]].insert(nearestSource[k][neighIdx]);
					}
				}
			}

			// Store in homogeneous datastructure
			int maxNeighNum = 0;
			for (int i = 0; i < neighborsList.size(); ++i) {
				if (neighborsList[i].size() > maxNeighNum) {
					maxNeighNum = neighborsList[i].size();
				}
			}
			neighLevelK.resize(DoF[k + 1], maxNeighNum);
			neighLevelK.setConstant(-1);
			for (int i = 0; i < neighborsList.size(); ++i) {
				neighLevelK(i, 0) = i;
				int iCounter = 1;
				for (int node : neighborsList[i]) {
					if (node == i)	continue;
					if (iCounter >= maxNeighNum) break;
					neighLevelK(i, iCounter) = node;
					iCounter++;
				}
			}
			hierarchyTiming["next_neighborhood"] += nextNeighTimer.get_elapsed_ms();

			plf::nanotimer nextPosTimer;
			nextPosTimer.start();
			if (verbose) std::cout << "Setting up the point locations for the next level\n";

			// Setting up the DoF for the next level
			// tempPoints are the centers of the voronoi cells, each row for each voronoi cells
			Eigen::MatrixXd tempPoints(DoF[k + 1], levelPoints.cols());
			tempPoints.setZero();
			if (nested) {
				for (int coarseIdx = 0; coarseIdx < DoF[k + 1]; ++coarseIdx) {
					tempPoints.row(coarseIdx) = levelPoints.row(samples[k][coarseIdx]);
				}
			} else {
				std::vector<int> clusterSizes(DoF[k + 1]);
				for (int fineIdx = 0; fineIdx < DoF[k]; ++fineIdx) {
					int coarseIdx = nearestSource[k][fineIdx];
					tempPoints.row(coarseIdx) = tempPoints.row(coarseIdx) + levelPoints.row(fineIdx);
					++clusterSizes[coarseIdx];
				}
				for (int coarseIdx = 0; coarseIdx < DoF[k + 1]; ++coarseIdx) {
					if (clusterSizes[coarseIdx] == 1) {
						tempPoints.row(coarseIdx) = levelPoints.row(samples[k][coarseIdx]);
						for (int neighIdx : neighborsList[coarseIdx]) {
							tempPoints.row(coarseIdx) = tempPoints.row(coarseIdx) + levelPoints.row(samples[k][neighIdx]);
						}
						tempPoints.row(coarseIdx) = tempPoints.row(coarseIdx) / (neighborsList[coarseIdx].size() + 1.);
					} else {
						tempPoints.row(coarseIdx) = tempPoints.row(coarseIdx) / clusterSizes[coarseIdx];
					}
				}
			}
			if (debug) levelV.push_back(tempPoints);
			hierarchyTiming["next_positions"] += nextPosTimer.get_elapsed_ms();

			plf::nanotimer triangleFindingTimer;
			triangleFindingTimer.start();

			// Create triangles for this level based on Voronoi cells
			std::vector<std::vector<int>> tris;
			tris.reserve(DoF[k + 1] * maxNeighNum);
			std::vector<std::vector<int>> connectedTris(DoF[k + 1]);
			std::vector<Eigen::RowVector3d> triNormals;
			triNormals.reserve(DoF[k + 1] * maxNeighNum);
			int triIdx = 0;
			for (int coarseIdx = 0; coarseIdx < DoF[k + 1]; ++coarseIdx) {
				// Iterate over delaunay triangles
				int v2Idx, v3Idx;
				for (auto it = neighborsList[coarseIdx].begin(); it != neighborsList[coarseIdx].end(); it++) {
					v2Idx = *it;
					// We iterate over the coarse indices in order,
					// so if the neighboring idx is lower then the current coarseIdx,
					// it must have been considered before and be part of a triangle.
					if (v2Idx < coarseIdx) continue;
					for (auto it2 = std::next(it); it2 != neighborsList[coarseIdx].end(); it2++) {
						v3Idx = *it2;
						if (v3Idx < coarseIdx) continue;
						if (neighborsList[v2Idx].find(v3Idx) != neighborsList[v2Idx].end()) {
							tris.push_back({coarseIdx, v2Idx, v3Idx});
							Eigen::RowVector3d e12 = tempPoints.row(v2Idx) - tempPoints.row(coarseIdx);
							Eigen::RowVector3d e13 = tempPoints.row(v3Idx) - tempPoints.row(coarseIdx);
							triNormals.push_back(e12.cross(e13).normalized());
							connectedTris[coarseIdx].push_back(triIdx);
							connectedTris[v2Idx].push_back(triIdx);
							connectedTris[v3Idx].push_back(triIdx);
							++triIdx;
						}
					}
				}
			}
			tris.shrink_to_fit();
			triNormals.shrink_to_fit();
			if (debug) allTriangles.push_back(tris);	
			hierarchyTiming["triangle_finding"] += triangleFindingTimer.get_elapsed_ms();

			plf::nanotimer triangleSelectionTimer;
			triangleSelectionTimer.start();

			// Create local triangulation on each cluster (centralized at sample i)
			int notrisfound = 0;
			int edgesfound = 0;
			int fallbackCount = 0;
			if (debug) noTriFoundMap.push_back(std::vector<int>(DoF[k]));

			for (int fineIdx = 0; fineIdx < DoF[k]; ++fineIdx) {
				Eigen::RowVector3d finePoint = levelPoints.row(fineIdx);
				int coarseIdx = nearestSource[k][fineIdx];
				Eigen::RowVector3d coarsePoint = tempPoints.row(coarseIdx);
				std::vector<double> weights;

				if (nested && samples[k][nearestSource[k][fineIdx]] == fineIdx) {
					AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, coarseIdx, 1.));
					continue;
				}

				if (neighborsList[coarseIdx].empty()) {
					// If the coarse point has no neighbors,
					// set the weight to 1 for the coarse point.
					// Note: this should not happen.
					AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, coarseIdx, 1.));
				} else if (neighborsList[coarseIdx].size() == 1) {
					// If the coarse point only has one neighbor, no triangle can be created.
					// Thus, the weights are distributed w.r.t. the distance to each coarse point.
					int neighIdx = *neighborsList[coarseIdx].begin();
					Eigen::RowVector3d neighPoint = tempPoints.row(neighIdx);

					// get the distance to the two neighboring centroids
					Eigen::RowVector3d e12 = neighPoint - coarsePoint;
					double e12Length = max(e12.norm(), 1e-8);
					double w2 = (levelPoints.row(fineIdx) - coarsePoint).dot(e12.normalized()) / e12Length;
					w2 = min(max(w2, 0.), 1.);
					double w1 = 1 - w2;

					switch (weightingScheme) {
						case BARYCENTRIC:
							AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, coarseIdx, w1));
							AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, neighIdx, w2));
							break;
						case UNIFORM:
							weights = uniformWeights(2);
							AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, coarseIdx, weights[0]));
							AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, neighIdx, weights[1]));
							break;
						case INVDIST:
							std::vector<int> endPoints = {coarseIdx, neighIdx};
							weights = inverseDistanceWeights(tempPoints, finePoint, endPoints);
							AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, coarseIdx, weights[0]));
							AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, neighIdx, weights[1]));
							break;
					}
				} else {
					// Only keep triangle with minimum distance
					double minDistToTriangle = std::numeric_limits<double>::max();
					Eigen::RowVector3d minBary = {0., 0., 0.};
					std::vector<int> minTri;
					bool triFound = false;

					// Values are positive if inside and negative if not
					// Float value represents distance
					std::map<int, float> insideEdge;

					// Iterate over all triangles
					for (int triIdx : connectedTris[coarseIdx]) {
						std::vector<int> tri = tris[triIdx];
						// Make sure that the coarseIdx is in position 0, while keeping orientation
						while(tri[0] != coarseIdx) std::rotate(tri.begin(), tri.begin() + 1, tri.end());

						Eigen::RowVector3d bary = {0., 0., 0.};
						// If the triangle contains the point, the distance is positive, else it's negative
						double distToTriangle = inTriangle(finePoint, tri, triNormals[triIdx], tempPoints, bary, insideEdge);	
						if (distToTriangle >= 0. && distToTriangle < minDistToTriangle) {
							triFound = true;
							minDistToTriangle = distToTriangle;
							minTri = tri;
							minBary = bary;
							break;
						}		
					}

					if (triFound) {
						switch (weightingScheme) {
							case BARYCENTRIC:
								AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, minTri[0], minBary(0)));
								AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, minTri[1], minBary(1)));
								AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, minTri[2], minBary(2)));
								break;
							case UNIFORM:
								weights = uniformWeights(3);
								AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, minTri[0], weights[0]));
								AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, minTri[1], weights[1]));
								AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, minTri[2], weights[2]));
								break;
							case INVDIST:
								weights = inverseDistanceWeights(tempPoints, finePoint, minTri);
								AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, minTri[0], weights[0]));
								AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, minTri[1], weights[1]));
								AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, minTri[2], weights[2]));
								break;
						}
					} else {
						bool edgeFound = false;
						double minEdge = std::numeric_limits<double>::max();
						int minEdgeIdx = 0;
						for (const auto& element : insideEdge) {
							const auto& key = element.first;
							const auto& value = element.second;
							if (value >= 0. && value < minEdge) {
								edgeFound = true;
								minEdge = value;
								minEdgeIdx = key;
								break;
							}
						}
						if (edgeFound) {
							++edgesfound;
							Eigen::RowVector3d p2 = tempPoints.row(minEdgeIdx);
							Eigen::RowVector3d e12 = p2 - coarsePoint;
							double e12Length = max(e12.norm(), 1e-8);
							double w2 = (finePoint - coarsePoint).dot(e12.normalized()) / e12Length;
							w2 = min(max(w2, 0.), 1.);
							double w1 = 1. - w2;

							switch (weightingScheme) {
								case BARYCENTRIC:
									AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, coarseIdx, w1));
									AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, minEdgeIdx, w2));
									break;
								case UNIFORM:
									weights = uniformWeights(2);
									AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, coarseIdx, weights[0]));
									AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, minEdgeIdx, weights[1]));
									break;
								case INVDIST:
									std::vector<int> endPointsEdge = {coarseIdx, minEdgeIdx};
									weights = inverseDistanceWeights(tempPoints, finePoint, endPointsEdge);
									AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, coarseIdx, weights[0]));
									AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, minEdgeIdx, weights[1]));
									break;
							}
						} else {
							// Use closest three
							std::vector<int> prolongFrom(3);
							prolongFrom[0] = coarseIdx;

							std::vector<VertexPair> pointsDistances;
							for (int j = 0; j < neighLevelK.cols(); ++j) {
								int neighIdx = neighLevelK(coarseIdx, j);
								if (neighIdx < 0 || neighIdx == coarseIdx) continue;
								VertexPair vp = {neighIdx, (finePoint - tempPoints.row(neighIdx)).norm()};
								pointsDistances.push_back(vp);
							}
							std::sort(pointsDistances.begin(), pointsDistances.end(), std::less<VertexPair>());
							for (int j = 1; j < 3; ++j) {
								prolongFrom[j] = pointsDistances[j - 1].vId;
							}
							std::vector<double> weights = inverseDistanceWeights(tempPoints, finePoint, prolongFrom);
							for (int j = 0; j < prolongFrom.size(); j++) {
								AllTriplet.push_back(Eigen::Triplet<double>(fineIdx, prolongFrom[j], weights[j]));
							}
							++fallbackCount;
						}
						++notrisfound;
					}
				}
			}
			if (verbose) cout << "Number of triangles unfound: " << notrisfound << endl;
			if (verbose) cout << "Number of edges found: " << edgesfound << endl;
			if (verbose) cout << "Percentage of fallback: " << (double) fallbackCount / (double) DoF[k] * 100 << endl;
			hierarchyTiming["triangle_selection"] += triangleSelectionTimer.get_elapsed_ms();

			levelPoints = tempPoints;

			Eigen::SparseMatrix<double> ULevel;
			ULevel.resize(DoF[k], DoF[k + 1]);
			ULevel.setFromTriplets(AllTriplet.begin(), AllTriplet.end());
			U.push_back(ULevel);
			AllTriplet.clear();
			AllTriplet.shrink_to_fit();
			++k;
		}
	}

	double MultigridSolver::inTriangle(const Eigen::RowVector3d& p, const std::vector<int>& tri, const Eigen::RowVector3d& triNormal, const Eigen::MatrixXd& pos, Eigen::RowVector3d& bary, std::map<int, float>& insideEdge)
	{
		Eigen::RowVector3d v1, v2, v3;
		v1 = pos.row(tri[0]);
		v2 = pos.row(tri[1]);
		v3 = pos.row(tri[2]);
		Eigen::RowVector3d v1ToP = p - v1;
		Eigen::RowVector3d e12 = v2 - v1;
		Eigen::RowVector3d e13 = v3 - v1;

		double distToTriangle = (p - v1).dot(triNormal);
		Eigen::RowVector3d pProjected = p - distToTriangle * triNormal;

		double doubleArea = (v2 - v1).cross(v3 - v1).dot(triNormal);
		bary(0) = (v3 - v2).cross(pProjected - v2).dot(triNormal) / doubleArea;
		bary(1) = (v1 - v3).cross(pProjected - v3).dot(triNormal) / doubleArea;
		bary(2) = 1. - bary(0) - bary(1);
		
		if (insideEdge.find(tri[1]) == insideEdge.end()) {
			insideEdge[tri[1]] = ((v1ToP) - ((v1ToP).dot(e12) * (e12))).norm();
		}
		if (insideEdge.find(tri[2]) == insideEdge.end()) {
			insideEdge[tri[2]] = ((v1ToP) - ((v1ToP).dot(e13) * (e13))).norm();
		}
		if (bary(0) < 0. || bary(1) < 0.) {
			insideEdge[tri[1]] = -1.;
		}
		if (bary(0) < 0. || bary(2) < 0.) {
			insideEdge[tri[2]] = -1.;
		}

		if (bary(0) >= 0. && bary(1) >= 0. && bary(2) >= 0.) {
			return abs(distToTriangle);
		}

		return -1.;
	}

	std::vector<double> MultigridSolver::uniformWeights(const int& n_points) {
		std::vector<double> weights(n_points);
		std::fill(weights.begin(), weights.end(), 1. / n_points);
		return weights;
	}

	std::vector<double> MultigridSolver::inverseDistanceWeights(const Eigen::MatrixXd& pos, const Eigen::RowVector3d& p, const std::vector<int>& edges) {
		double sumWeight = 0.;
		std::vector<double> weights(edges.size());
		for (int j = 0; j < edges.size(); ++j) {
			weights[j] = 1. / max(1e-8, (p - pos.row(edges[j])).norm());
			sumWeight += weights[j];
		}
		for (int j = 0; j < weights.size(); ++j) {
			weights[j] = weights[j] / sumWeight;
		}
		return weights;
	}

	double MultigridSolver::computeAverageEdgeLength(const Eigen::MatrixXd& pos, const Eigen::MatrixXi& neigh) {
		double sumLength = 0;
		int nEdges = 0;
		for (int i = 0; i < pos.rows(); ++i) {
			Eigen::Vector3d p1 = pos.row(i);
			for (int j = 0; j < neigh.cols(); ++j) {
				if (neigh(i, j) < 0) continue;
				Eigen::Vector3d p2 = pos.row(neigh(i, j));
				double dist = (p1 - p2).norm();
				if (dist > 0) {
					sumLength += dist;
					++nEdges;
				}
			}
		}
		return sumLength / (double) nEdges;
	}

	std::vector<int> MultigridSolver::maximumDeltaIndependentSet(const Eigen::MatrixXd& pos, const Eigen::MatrixXi& edges, const double& radius) {
		std::vector<bool> visited(edges.rows());
		std::vector<int> selection;
		for (int i = 0; i < edges.rows(); ++i) {
			if (!visited[i]) {
				selection.push_back(i);
				for (int j = 0; j < edges.cols(); ++j) {
					int neighIdx = edges(i, j);
					if (neighIdx < 0) break;
					double dist = (pos.row(i) - pos.row(neighIdx)).norm();
					if (dist < radius) {
						visited[neighIdx] = true;
					}
				}
			}
		}
		return selection;
	}

	std::vector<int> MultigridSolver::maximumDeltaIndependentSet(const Eigen::MatrixXd& pos, const Eigen::MatrixXi& edges, const double& radius, Eigen::VectorXd& D, std::vector<size_t>& nearestSourceK) {
		std::vector<bool> visited(edges.rows());
		std::vector<int> selection;
		int sampleIdx = 0;
		for (int i = 0; i < edges.rows(); ++i) {
			if (!visited[i]) {
				selection.push_back(i);
				nearestSourceK[i] = sampleIdx;
				for (int j = 0; j < edges.cols(); ++j) {
					int neighIdx = edges(i, j);
					if (neighIdx < 0) break;
					double dist = (pos.row(i) - pos.row(neighIdx)).norm();
					if (dist < radius) {
						visited[neighIdx] = true;
						if (dist < D(neighIdx)) {
							D(neighIdx) = dist;
							nearestSourceK[neighIdx] = sampleIdx;
						}
					}
				}
				++sampleIdx;
			}
		}
		return selection;
	}

	std::vector<int> MultigridSolver::fastDiskSample(const Eigen::MatrixXd& pos, const Eigen::MatrixXi& edges, const double& radius, Eigen::VectorXd& D, std::vector<size_t>& nearestSourceK) {
		std::vector<bool> visited(edges.rows());
		std::vector<int> selection;
		int sampleIdx = 0;
		for (int i = 0; i < edges.rows(); ++i) {
			if (!visited[i]) {
				selection.push_back(i);
				nearestSourceK[i] = sampleIdx;
				for (int j = 0; j < edges.cols(); ++j) {
					int neighIdx = edges(i, j);
					if (neighIdx < 0) break;
					double dist = (pos.row(i) - pos.row(neighIdx)).norm();
					if (dist < radius) {
						visited[neighIdx] = true;
						if (dist < D(neighIdx)) {
							D(neighIdx) = dist;
							nearestSourceK[neighIdx] = sampleIdx;
						}
						for (int j2 = 0; j2 < edges.cols(); ++j2) {
							int neighIdx2 = edges(neighIdx, j2);
							if (neighIdx2 < 0) break;
							double dist2 = dist + (pos.row(neighIdx) - pos.row(neighIdx2)).norm();
							if (dist2 < radius) {
								visited[neighIdx2] = true;
								if (dist2 < D(neighIdx2)) {
									D(neighIdx2) = dist2;
									nearestSourceK[neighIdx2] = sampleIdx;
								}
							}
						}
					}
				}
				++sampleIdx;
			}
		}
		return selection;
	}

	void MultigridSolver::constructDijkstraWithCluster(const Eigen::MatrixXd& points, const std::vector<int>& source, const Eigen::MatrixXi& neigh, int k, Eigen::VectorXd& D, std::vector<size_t>& nearestSourceK)
	{
		std::priority_queue<VertexPair, std::vector<VertexPair>, std::greater<VertexPair>> DistanceQueue;
		if (nearestSourceK.empty())	nearestSourceK.resize(points.rows(), source[0]);

		for (int i = 0; i < source.size(); ++i) {
			D(source[i]) = 0.0;
			VertexPair vp{ source[i], D(source[i]) };
			DistanceQueue.push(vp);
			nearestSourceK[source[i]] = i;
		}

		int curSource;
		while (!DistanceQueue.empty()) {
			VertexPair vp1 = DistanceQueue.top();
			curSource = nearestSourceK[vp1.vId];
			Eigen::RowVector3d vertex1 = points.row(vp1.vId);
			DistanceQueue.pop();

			//for (int vNeigh : neigh.row(vp1.vId)) {
			for (int i = 0; i < neigh.cols(); ++i) {
				int vNeigh = neigh(vp1.vId, i);

				if (vNeigh >= 0) {
					double dist, distTemp;
					Eigen::RowVector3d vertex2 = points.row(vNeigh);
					dist = (vertex2 - vertex1).norm();
					distTemp = vp1.distance + dist;
					if (distTemp < D(vNeigh)) {
						// Assign a new distance
						D(vNeigh) = distTemp;
						VertexPair v2{ vNeigh, distTemp };
						DistanceQueue.push(v2);


						// Assign the nearest source to a certain point
						nearestSourceK[vNeigh] = curSource;
					}
				}
			}
		}
	}


	double  MultigridSolver::multiGridVCycleGS(Eigen::SparseMatrix<double>& A, Eigen::MatrixXd& b, Eigen::MatrixXd& x, int k, bool isDebug)
	{
		double tolPre = 1e-15, tolPost = 1e-15;
		//[1] Smoothing
		GaussSeidelSmoother(A, b, x, preIters, tolPre, isDebug);

		//[2] Residual
		Eigen::MatrixXd res = b - A * x;

		//[3] Restriction
		Eigen::MatrixXd resRest = U[k].transpose() * res;

		//[4] Recursion
		Eigen::MatrixXd eps;
		eps.setZero(resRest.rows(), resRest.cols());
		if (k == U.size() - 1) {
			eps = coarsestSolver.solve(resRest);
		}
		else {
			multiGridVCycleGS(Abar[k + 1], resRest, eps, k + 1, isDebug);
		}

		//[5] Prolongation
		x = x + U[k] * eps;

		//[6] Post-Smoothing
		GaussSeidelSmoother(A, b, x, postIters, tolPost, isDebug);

		return 0.0;
	}

	//!< Perform Multigrid using F-Cycle and Gauss-Seidel smoothing
	double  MultigridSolver::multiGridFCycleGS(Eigen::SparseMatrix<double>& A, Eigen::MatrixXd& b, Eigen::MatrixXd& x, int k, bool isDebug)
	{
		double tolPre = 1e-15, tolPost = 1e-15;
		//[1] Smoothing
		GaussSeidelSmoother(A, b, x, preIters, tolPre, isDebug);

		//[2] Residual
		Eigen::MatrixXd res = b - A * x;

		//[3] Restriction
		Eigen::MatrixXd resRest = U[k].transpose() * res;

		//[4] Recursion
		Eigen::MatrixXd eps;
		eps.setZero(resRest.rows(), resRest.cols());
		if (k == U.size() - 1) {
			eps = coarsestSolver.solve(resRest);
		}
		else {
			multiGridFCycleGS(Abar[k + 1], resRest, eps, k + 1, isDebug);
		}

		//[5] Prolongation
		x = x + U[k] * eps;

		//[6] Post-Smoothing
		GaussSeidelSmoother(A, b, x, postIters, tolPost, isDebug);

		//[7] Residual
		res = b - A * x;

		//[8] Restriction
		resRest = U[k].transpose() * res;

		//[9] Recursion
		if (k == DoF.size() - 2) {
			eps = coarsestSolver.solve(resRest);
		}
		else {
			multiGridVCycleGS(Abar[k + 1], resRest, eps, k + 1, isDebug);
		}

		//[10] Prolongation
		x = x + U[k] * eps;

		//[11] Post-Smoothing
		GaussSeidelSmoother(A, b, x, postIters, tolPost, isDebug);

		return 0.0;
	}

	//!< Perform Multigrid using W-Cycle and Gauss-Seidel smoothing
	double  MultigridSolver::multiGridWCycleGS(Eigen::SparseMatrix<double>& A, Eigen::MatrixXd& b, Eigen::MatrixXd& x, int k, bool isDebug)
	{
		double tolPre = 1e-15, tolPost = 1e-15;
		//[1] Smoothing
		GaussSeidelSmoother(A, b, x, preIters, tolPre, isDebug);

		//[2] Residual
		Eigen::MatrixXd res = b - A * x;

		//[3] Restriction
		Eigen::MatrixXd resRest = U[k].transpose() * res;

		//[4] Recursion
		Eigen::MatrixXd eps;
		eps.setZero(resRest.rows(), resRest.cols());
		if (k == U.size() - 1) {
			eps = coarsestSolver.solve(resRest);
		}
		else {
			multiGridWCycleGS(Abar[k + 1], resRest, eps, k + 1, isDebug);
		}

		//[5] Prolongation
		x = x + U[k] * eps;

		//[6] Post-Smoothing
		GaussSeidelSmoother(A, b, x, postIters, tolPost, isDebug);

		//[7] Residual
		res = b - A * x;

		//[8] Restriction
		resRest = U[k].transpose() * res;

		//[9] Recursion
		if (k == DoF.size() - 2) {
			eps = coarsestSolver.solve(resRest);
		}
		else {
			multiGridWCycleGS(Abar[k + 1], resRest, eps, k + 1, isDebug);
		}

		//[10] Prolongation
		x = x + U[k] * eps;

		//[11] Post-Smoothing
		GaussSeidelSmoother(A, b, x, postIters, tolPost, isDebug);

		return 0.0;
	}

	void MultigridSolver::GaussSeidelSmoother(Eigen::SparseMatrix<double>& LHS, Eigen::MatrixXd& rhs,
		Eigen::MatrixXd& x, int maxIter, double tol, bool isDebug)
	{
		int dim = x.cols();
		if (dim == 1) {
			for (int i = 0; i < maxIter; i++) {
				for (int k = 0; k < LHS.outerSize(); ++k) {
					double sum = 0.0;
					for (Eigen::SparseMatrix<double>::InnerIterator it(LHS, k); it; ++it) {
						if (it.row() != k) {
							sum += it.value() * x(it.row());
						}
					}
					x(k) = (rhs(k) - sum) / LHS.coeffRef(k, k);
				}
			}
		}
		else {
			for (int i = 0; i < maxIter; i++) {
				for (int iCol = 0; iCol < dim; ++iCol) {
					for (int k = 0; k < LHS.outerSize(); ++k) {
						double sum = 0.0;
						for (Eigen::SparseMatrix<double>::InnerIterator it(LHS, k); it; ++it) {
							if (it.row() != k) {
								sum += it.value() * x(it.row(), iCol);
							}
						}
						x(k, iCol) = (rhs(k, iCol) - sum) / LHS.coeffRef(k, k);
					}
				}
			}
		}
	}

	double MultigridSolver::residualCheck(const Eigen::SparseMatrix<double>& A, const Eigen::MatrixXd& b, const Eigen::MatrixXd& x, int type)
	{
		// 0: rel. norm (||Ax-b||/||b||)		1: L2 M^-1 (||Ax-b||_{M-1}/||b||_{M-1})		2: L2 M (||Ax-b||{M}/||b||_{M})		3: Abs (Ax-b).norm()
		double residue;
		switch (type)
		{
		case 0:
		{
			Eigen::VectorXd absResidue(b.cols());
			for (int i = 0; i < b.cols(); ++i) {
				absResidue(i) = (A * x.col(i) - b.col(i)).norm() / b.col(i).norm();
			}
			residue = absResidue.maxCoeff();
			break;
		}
		case 1:
		{
			Eigen::VectorXd absResidue(b.cols());
			Eigen::VectorXd residual;
			for (int i = 0; i < b.cols(); ++i) {
				residual = A * x.col(i) - b.col(i);
				double n1 = residual.transpose() * Minv * residual;
				double n2 = b.col(i).transpose() * Minv * b.col(i);
				absResidue(i) = sqrt(n1 / n2);
			}
			residue = absResidue.maxCoeff();
			break;
		}
		case 2:
		{
			Eigen::VectorXd absResidue(b.cols());
			Eigen::VectorXd residual;
			for (int i = 0; i < b.cols(); ++i) {
				residual = A * x.col(i) - b.col(i);
				double n1 = residual.transpose() * M * residual;
				double n2 = b.col(i).transpose() * M * b.col(i);
				absResidue(i) = sqrt(n1 / n2);
			}
			residue = absResidue.maxCoeff();
			break;
		}
		case 3:
			residue = (A * x - b).norm();
			break;
		default:
			break;
		}

		return residue;
	}

	void MultigridSolver::solve(Eigen::SparseMatrix<double>& LHS, Eigen::MatrixXd& rhs, Eigen::MatrixXd& x)
	{
		std::chrono::high_resolution_clock::time_point t0, t1, t2;
		std::chrono::duration<double> duration1, duration2;
		convergence.reserve(maxIter);

		Eigen::MatrixXd mx = rhs;

		if (verbose) std::cout << "Solve the linear system using Gravo MG solver \n";

		double residue = std::numeric_limits<double>::max();
		Eigen::VectorXd rhs0;
		Eigen::VectorXd rhs1;
		Eigen::VectorXd rhs2;
		rhs0 = mx.col(0);
		if (mx.cols() > 1) {
			rhs1 = mx.col(1);
			rhs2 = mx.col(2);
		}

		if (isSmootherGaussSeidel)
		{
			plf::nanotimer solverTotalTimer, reductionTimer;
			solverTotalTimer.start(); reductionTimer.start();

			if (verbose) cout << "Reducing system" << endl;

			Abar.resize(U.size() + 1);
			//Mbar.resize(DoF.size());
			Abar[1] = U[0].transpose() * LHS * U[0];
			for (int k = 2; k < U.size() + 1; ++k) {
				Abar[k] = U[k - 1].transpose() * Abar[k - 1] * U[k - 1];
			}

			solverTiming["reduction"] = reductionTimer.get_elapsed_ms();

			plf::nanotimer coarseSolveTimer;
			coarseSolveTimer.start();

			if (verbose) cout << "Solving reduced system" << endl;

			coarsestSolver.compute(Abar[U.size()]);

			solverTiming["coarsest_solve"] = coarseSolveTimer.get_elapsed_ms();

			plf::nanotimer cyclesTimer;
			cyclesTimer.start();
			int iterCount = 0;
			if (cycleType == 0) {
				if (verbose) std::cout << "V-CYCLE \n";
				Eigen::Vector3d relResidue;
				do {
					residue = multiGridVCycleGS(LHS, mx, x, 0, false);
					residue = residualCheck(LHS, mx, x, stoppingCriteria);
					convergence.push_back({cyclesTimer.get_elapsed_ms(), residue});
					if (verbose) printf("%d,%f,%.14f \n", ++iterCount, cyclesTimer.get_elapsed_ms(), residue);
					else ++iterCount;
				} while ((residue > accuracy) && (iterCount < maxIter));
				convergence.shrink_to_fit();
			}
			else if (cycleType == 1) {
				if (verbose) std::cout << "F-CYCLE \n";
				Eigen::Vector3d relResidue;
				do {
					residue = multiGridFCycleGS(LHS, mx, x, 0, false);
					residue = residualCheck(LHS, mx, x, stoppingCriteria);
					if (verbose) printf("Iteration %d: residue=%.10f \n", ++iterCount, residue);
					else ++iterCount;
				} while ((residue > accuracy) && (iterCount < maxIter));
			}
			else if (cycleType == 2) {
				if (verbose) std::cout << "W-CYCLE \n";
				Eigen::Vector3d relResidue;
				do {
					residue = multiGridWCycleGS(LHS, mx, x, 0, false);
					residue = residualCheck(LHS, mx, x, stoppingCriteria);
					if (verbose) printf("Iteration %d: residue=%.10f \n", ++iterCount, residue);
					else ++iterCount;
				} while ((residue > accuracy) && (iterCount < maxIter));
			}
			else {
				if (verbose) std::cout << "ERROR! The cycle type should only be 0-2! \n";
				return;
			}

			solverTiming["cycles"] = cyclesTimer.get_elapsed_ms();
			solverTiming["solver_total"] = solverTotalTimer.get_elapsed_ms();
			solverTiming["iterations"] = iterCount * 1.0;
			solverTiming["residue"] = residue;
		}
			
	}

}

