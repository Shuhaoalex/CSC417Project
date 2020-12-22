#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <cmath>
#include "Grid3D.h"
#include "utils.h"
#include "constants.h"
#include "Particle.h"

Grid3D grid(grid_size, h);
std::map<std::pair<size_t, size_t>, Spring> springs;
std::vector<Particle> particle_list;
constants c;
bool animating = false;

void reset_c() {
    c.g = Eigen::Vector3d(0,-0.005,0);
    c.num_particles = 10000;

    c.dt = 1;

    // double density paramter
    // c.k = 0.0004;
    // c.k_near = 0.001;
    c.k = 0.001;
    c.k_near = 0.001;
    c.rho0 = 10;

    // viscosity parameter
    c.sigma = 0.1;
    c.beta = 1;


    c.collision_e = 0.3;

    // plasticity parameter
    c.alpha = 0.3;
    c.gamma = 0.1;
    c.k_spring = 0.3;
}

void sphere_setup() {
    reset_c();
    animating = false;
    grid.clear();
    springs.clear();
    particle_list.clear();
    particle_list.reserve(c.num_particles);
    for (int i = 0; i < c.num_particles; ++i) {
        Eigen::Vector3d currx = Eigen::Vector3d::Random();
        while (currx.norm() > 1) {
            currx = Eigen::Vector3d::Random();
        }
        currx = (currx + Eigen::Vector3d(1,1,1))*box_size/2.;
        particle_list.emplace_back(currx, Eigen::Vector3d(0,0,0), i);
    }
    grid.insert_particles(particle_list);
}

void add_sphere() {
    int num_particles = 500;
    particle_list.reserve(particle_list.size() + num_particles);
    for (int i = 0; i < num_particles; ++i) {
        Eigen::Vector3d currx = Eigen::Vector3d::Random();
        while (currx.norm() > 1) {
            currx = Eigen::Vector3d::Random();
        }
        currx *= 0.35;
        currx = (currx + Eigen::Vector3d(1,1.5,1))*box_size/2.;
        particle_list.emplace_back(currx, Eigen::Vector3d(0,0,0), i);
    }
    grid.clear();
    grid.insert_particles(particle_list);
}

void weightless_box_setup() {
    reset_c();
    animating = false;
    c.k = 0.002;
    c.k_near = 0.002;
    grid.clear();
    springs.clear();
    particle_list.clear();
    particle_list.reserve(c.num_particles);
    c.g = Eigen::Vector3d::Zero();
    for (int i = 0; i < c.num_particles; ++i) {
        Eigen::Vector3d currx = Eigen::Vector3d::Random() / 2 * 0.8;
        currx = (currx + Eigen::Vector3d(0.5,0.5,0.5)) * box_size;
        particle_list.emplace_back(currx, Eigen::Vector3d(0,0,0), i);
    }
    grid.insert_particles(particle_list);
}

void weightless_box_setup2() {
    reset_c();
    animating = false;
    c.k = 0.0008;
    c.k_near = 0.001;
    c.num_particles = 500;
    grid.clear();
    springs.clear();
    particle_list.clear();
    particle_list.reserve(c.num_particles);
    c.g = Eigen::Vector3d::Zero();
    for (int i = 0; i < c.num_particles; ++i) {
        Eigen::Vector3d currx = Eigen::Vector3d::Random() / 2;
        currx.tail<2>() *= 0.3;
        currx = (currx + Eigen::Vector3d(0.5,0.5,0.5)) * box_size;
        particle_list.emplace_back(currx, Eigen::Vector3d(0,0,0), i);
    }
    grid.insert_particles(particle_list);
}

bool key_down_callback(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifiers) {
    switch(key) {
        case 'A':
            animating = !animating;
            break;
        case 'Q':
            sphere_setup();
            break;
        case 'W':
            weightless_box_setup();
            break;
        case 'E':
            weightless_box_setup2();
            break;
        case 'S':
            add_sphere();
            break;
    }
    return false;
}


int main(int argc, char** argv) {
    using namespace Eigen;
    igl::opengl::glfw::Viewer v;
    v.callback_key_down = key_down_callback;
    RowVector3d color(103, 51, 32);
    color /= 255;
    weightless_box_setup2();
    // while (true) {
    //     simulation_step(particle_list, grid);
    // }

    // checkout this tutorial
    // https://stackoverflow.com/questions/55878584/clear-move-animate-a-point-on-libigl-viewer
    v.callback_pre_draw = [&](igl::opengl::glfw::Viewer & )->bool
    {
        // before running next time step, clear previous pts.
        v.data().clear_points();
        if (animating) simulation_step(particle_list, grid, springs, c);
        // add points one by one
        for (auto& par:particle_list) {
            v.data().add_points(par.x.transpose().array() - 0.5 * box_size, color);
        }
        return false;
    };

    v.data().point_size = 4;
    v.core().is_animating = true;
    v.launch();
    return 0;
}

