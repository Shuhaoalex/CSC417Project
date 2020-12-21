#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include "Grid3D.h"
#include "utils.h"
#include "constants.h"
#include "Particle.h"

Grid3D grid(grid_size, h);
std::vector<Particle> particle_list;

void setup() {
    std::cout << "start setting up\n";
    particle_list.clear();
    particle_list.reserve(num_particles);
    for (int i = 0; i < num_particles; ++i) {
        Eigen::Vector3d currx = Eigen::Vector3d::Random();
        currx = (currx + Eigen::Vector3d(1,1,1))*box_size/2.;
        particle_list.emplace_back(currx, Eigen::Vector3d(0,0,0), i);
    }
    grid.insert_particles(particle_list);
    std::cout << "finished setting up\n";
}

int main(int argc, char** argv) {
    using namespace Eigen;
    igl::opengl::glfw::Viewer v;
    RowVector3d color(103, 51, 32);
    color /= 255;
    setup();
    // while (true) {
    //     simulation_step(particle_list, grid);
    // }

    // checkout this tutorial
    // https://stackoverflow.com/questions/55878584/clear-move-animate-a-point-on-libigl-viewer
    v.callback_pre_draw = [&](igl::opengl::glfw::Viewer & )->bool
    {
        // before running next time step, clear previous pts.
        v.data().clear_points();

        simulation_step(particle_list, grid);
        // add points one by one
        for (auto& par:particle_list) {
            v.data().add_points(par.x.transpose(), color);
        }
        return false;
    };

    v.data().point_size = 4;
    v.core().is_animating = true;
    v.launch();
    return 0;
}

