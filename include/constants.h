#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__
#include <Eigen/Core>
#include <cmath>
const double h = 0.1;
const int grid_size = 10;
const double box_size = h * grid_size;
const int num_particles = 10000;
const double k = 0.000005;
const double dt = 1;
const double k_near = k * 100;
const double rho0 = 3.0;
const double sigma = 0.1;
const double beta = 1;
const double collision_e = 1.0;
// const double alpha = 1.0;
// const double gamma = 0.1;
const Eigen::Vector3d g = Eigen::Vector3d(0,-0.005,0);
#endif