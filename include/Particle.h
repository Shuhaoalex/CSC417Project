#ifndef __PARTICLE_H__
#define __PARTICLE_H__
#include <Eigen/Core>
#include <vector>

struct Particle {
    size_t id;
    Eigen::Vector3d x_prev;
    Eigen::Vector3d x;
    Eigen::Vector3d v;
    Eigen::Vector3i grid_loc;
    std::vector<Particle*> neighbors;
    inline Particle(Eigen::Vector3d _x, Eigen::Vector3d _v, size_t _id):x(_x), v(_v), id(_id){};
};

#endif