#include "gravomg/utility.h"

#include <igl/doublearea.h>

#include <fstream>


namespace MGBS {
	void scaleMesh(Eigen::MatrixXd& V, const Eigen::MatrixXi& F, double scaleRatio)
	{
		//cout << "> Scaling the mesh to have a unit length \n";
		Eigen::Vector3d minV;
		Eigen::Vector3d maxV;
		Eigen::Vector3d length;
		Eigen::MatrixXd MatSubs(V.rows(), V.cols());
		Eigen::MatrixXd MatAdd(V.rows(), V.cols());
		double maxVal;
		double scalingRatio = scaleRatio; //dia-meter, not radius


		/* Get the min and max coefficients */
		for (int i = 0; i < V.cols(); i++)
		{
			minV(i) = V.col(i).minCoeff();
			maxV(i) = V.col(i).maxCoeff();
			length(i) = maxV(i) - minV(i);
			MatSubs.col(i).setConstant(minV(i));
		}

		maxVal = length.maxCoeff();

		/* Translate to the Origin */
		V = V - MatSubs;

		/* Scale w.r.t the longest axis */
		V = V * (scalingRatio / maxVal);

		for (int i = 0; i < V.cols(); i++) {
			maxV(i) = V.col(i).maxCoeff();
			MatAdd.col(i).setConstant(0.5 * maxV(i));
		}

		/* Translate s.t. the center is in the Origin */
		V = V - MatAdd;

		for (int i = 0; i < V.cols(); i++)
		{
			minV(i) = V.col(i).minCoeff();
			maxV(i) = V.col(i).maxCoeff();
			length(i) = maxV(i) - minV(i);
		}
		maxVal = length.maxCoeff();
	}

	void normalize_unit_area(Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
	{
		Eigen::VectorXd FA;
		igl::doublearea(V, F, FA);

		double scale = sqrt(FA.sum() / 2);
		V /= scale;
				
		V.col(0) = V.col(0).array() - V.col(0).mean();
		V.col(1) = V.col(1).array() - V.col(1).mean();
		//V.col(2) = V.col(2).array() - V.col(2).minCoeff();
		V.col(2) = V.col(2).array() - V.col(2).mean();

		std::cout << "Check mass matrix after normalization \n";
		igl::doublearea(V, F, FA);
		std::cout << "Area: " << FA.sum() / 2 << " and M(0,0)=" << FA(0) / 2 << std::endl;
	}

}
