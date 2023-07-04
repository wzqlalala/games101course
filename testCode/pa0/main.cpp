#include<cmath>
#include"Eigen3/Eigen/Core"
#include"Eigen3/Eigen/Dense"
#include "Eigen3/Eigen/Eigen"
#include<iostream>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace Eigen;
int main(){

    // Basic Example of cpp
    std::cout << "Example of cpp \n";
    float a = 1.0, b = 2.0;
    std::cout << a << std::endl;
    std::cout << a/b << std::endl;
    std::cout << std::sqrt(b) << std::endl;
    std::cout << std::acos(-1) << std::endl;
    std::cout << std::sin(30.0/180.0*acos(-1)) << std::endl;

    // Example of vector
    std::cout << "Example of vector \n";
    // vector definition
    Eigen::Vector3f v(1.0f,2.0f,3.0f);
    Eigen::Vector3f w(1.0f,0.0f,0.0f);
    // vector output
    std::cout << "Example of output \n";
    std::cout << v << std::endl;
    // vector add
    std::cout << "Example of add \n";
    std::cout << v + w << std::endl;
    // vector scalar multiply
    std::cout << "Example of scalar multiply \n";
    std::cout << v * 3.0f << std::endl;
    std::cout << 2.0f * v << std::endl;

    // Example of matrix
    std::cout << "Example of matrix \n";
    // matrix definition
    Eigen::Matrix3f i,j;
    i << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
    j << 2.0, 3.0, 1.0, 4.0, 6.0, 5.0, 9.0, 7.0, 8.0;
    // matrix output
    std::cout << "Example of output \n";
    std::cout << i << std::endl;
    // matrix add i + j
    // matrix scalar multiply i * 2.0
    // matrix multiply i * j
    // matrix multiply vector i * v

	Vector3f vec2(2, 1, 1); float theta = 45 * M_PI / 180.0;
	Matrix3f model1;
	model1 << cos(theta), sin(theta), 1, -sin(theta), cos(theta), 2, 0, 0, 1;
	Vector3f res = model1 * vec2;
    return 0;
}