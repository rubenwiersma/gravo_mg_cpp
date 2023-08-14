#ifndef UTILITY_H
#define UTILITY_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "plf_nanotimer.h"

#include <set>
#include <vector>
#include <random>

#include <Eigen/Eigen>

using namespace std;

namespace GravoMG {
	void scaleMesh(Eigen::MatrixXd& V, const Eigen::MatrixXi& F, double scaleRatio = 1.0);
}

#endif // !UTILITY_H