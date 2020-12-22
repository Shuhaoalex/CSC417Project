#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__
#include <Eigen/Core>

const double h = 0.1;
const int grid_size = 10;
const double box_size = h * grid_size;
struct constants {
    
    Eigen::Vector3d g = Eigen::Vector3d(0,-0.005,0);
    int num_particles = 10000;

    double dt = 1;

    // double density paramter
    double k = 0.00001;
    double k_near = 0.005;
    double rho0 = 10;

    // viscosity parameter
    double sigma = 0.1;
    double beta = 1;


    double collision_e = 0.3;

    // plasticity parameter
    double alpha = 0.3;
    double gamma = 0.1;
    double k_spring = 0.3;
};
// const double L = h;

#endif