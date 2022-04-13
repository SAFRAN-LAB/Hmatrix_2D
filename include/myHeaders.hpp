#ifndef myHeaders_HPP
#define myHeaders_HPP

#include <iostream>
#include <chrono>
#include <ctime>
#include <set>
#include <string>
#include <vector>
#include <random>
#include <cmath>
#include <complex>
#include <iterator>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <list>
#include <algorithm>
#include <numeric>
#include <stack>

#define _USE_MATH_DEFINES
#include <Eigen/Dense>
#include <Eigen/Sparse>

using std::string;
using namespace Eigen;

#ifdef USE_DOUBLE
using dtype=double;
using dtype_base=double;
using Mat=Eigen::MatrixXd;
using Vec=Eigen::VectorXd;
#define abs_(x) ((x<0) ? (-x) : (x))
#define conj_(x) (x)
#define Calc_dist(x1,y1,x2,y2) double(sqrt(pow(x1-x2,2)+pow(y1-y2,2)))
#endif

#ifdef USE_COMPLEX64
using dtype=std::complex<double>;
using dtype_base=double;
using Mat=Eigen::MatrixXcd;
using Vec=Eigen::VectorXcd;
#define abs_(x) (std::abs((x)))
#define conj_(x) (std::conj((x)))
const std::complex<double> I(0.0, 1.0);
#endif


#ifdef USE_FLOAT
using dtype=float;
using dtype_base=float;
using Mat=Eigen::MatrixXf;
using Vec=Eigen::VectorXf;
#define Calc_dist(x1,y1,x2,y2) float(sqrt(pow(x1-x2,2)+pow(y1-y2,2)))
#endif

#ifdef USE_COMPLEX32
using dtype=std::complex<float>;
using dtype_base=float;
using Mat=Eigen::MatrixXcf;
using Vec=Eigen::VectorXcf;
const std::complex<float> I(0.0, 1.0);
#endif

std::string inline getTimeStamp()
{
        std::string s;
        s = "\n" + string(__DATE__) + "," + string(__TIME__) + "\n";
        return s;
}

namespace Vec_ops
{
dtype_base inline relative_error(Vec& X_ori,Vec& X_comp)
{
        dtype_base result = 0.0;
        result = (X_ori - X_comp).norm() / X_ori.norm();
        return result;
}
}

#endif
//const std::complex<double> I(0.0, 1.0);
