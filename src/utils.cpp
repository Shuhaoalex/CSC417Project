#include "utils.h"
#include "constants.h"
#include <iostream>

void retrieve_neiboroughs(std::vector<Particle> &pl, Grid3D & grid) {
    grid.update_particles_pos();
    #pragma omp parallel for
    for (int i = 0; i < pl.size(); ++i) {
        grid.get_neighbors(pl[i], pl[i].neighbors);
    }
}

void double_density_relaxation(std::vector<Particle> &pl) {
    #pragma omp parallel for
    for (int i = 0; i < pl.size(); ++i) {
        auto& par = pl[i];
        double rho = 0;
        double rho_near = 0;
        for (auto& np: par.neighbors) {
            Eigen::Vector3d r = np->x - par.x;
            double q = r.norm() / h;
            if (q < 1) {
                rho += pow(1-q, 2);
                rho_near += pow(1-q, 3);
            }
        }
        double P = k * (rho - rho0);
        double P_near = k_near * rho_near;
        Eigen::Vector3d dx = Eigen::Vector3d::Zero();
        for (auto& np: par.neighbors) {
            Eigen::Vector3d r = np->x - par.x;
            double q = r.norm() / h;
            if (q < 1) {
                Eigen::Vector3d D = pow(dt, 2) * (P * (1-q) + P_near * pow(1-q, 2)) / 2 * r.normalized();
                for (int j = 0; j < 3; ++j) {
                    #pragma omp atomic update
                    np->x(j) += D(j);
                }
                dx -= D;
            }
        }
        par.x += dx;
    }
}

void viscosity_impulses(std::vector<Particle> &pl) {
    #pragma omp parallel for
    for (int i = 0; i < pl.size(); ++i) {
        auto& par = pl[i];
        for (auto& np:par.neighbors) {
            if (par.id < np->id) {
                Eigen::Vector3d r = np->x - par.x;
                double q = r.norm() / h;
                if (q < 1) {
                    double u = (par.v - np->v).dot(r.normalized());
                    if (u > 0) {
                        Eigen::Vector3d I = dt * (1-q) * (sigma * u + beta * u * u) / 2 * r.normalized();
                        for (int j = 0; j < 3; ++j) {
                            #pragma omp atomic update
                            par.v(j) -= I(j);
                            #pragma omp atomic update
                            np->v(j) += I(j);
                        }
                    }
                }
            }
        }
    }
}

void apply_gravity(std::vector<Particle> &pl) {
    #pragma omp parallel for
    for (int i = 0; i < pl.size(); ++i) {
        pl[i].v += dt * g;
    }
}

void avection(std::vector<Particle> &pl) {
    #pragma omp parallel for
    for (int i = 0; i < pl.size(); ++i) {
        pl[i].x_prev = pl[i].x;
        pl[i].x += dt * pl[i].v;
    }
}

void update_vs(std::vector<Particle> &pl) {
    #pragma omp parallel for
    for (int i = 0; i < pl.size(); ++i) {
        auto&p = pl[i];
        p.v = (p.x - p.x_prev) / dt;
        for (int ax = 0; ax < 3; ++ax) {
            if (p.x(ax) < 0) {
                p.x(ax) = 0;
                p.v(ax) *= -collision_e;
            } else if (p.x(ax) > box_size) {
                p.x(ax) = box_size;
                p.v(ax) *= -collision_e;
            }
        }
    }
}

// void adjustSprings(std::vector<Particle> &pl, std::vector<std::pair<size_t, size_t>> &springs) {
//     for (auto& par:pl) {
//         for (auto& np:par.neighbors) {
//             if (par.id < np->id) {
//                 Eigen::Vector3d r = np->x - par.x;
//                 double q = r.norm() / h;
                
//             }
//         }
//     }
// }

void adjustSprings(std::vector<Particle> &pl, std::map<std::pair<size_t, size_t>, double> &springs) {
    for (size_t i = 0; i < pl.size(); ++i) {
        auto& par = pl[i];
        for (auto& np: par.neighbors) {
            Eigen::Vector3d r = np->x - par.x;
            double 
            double q = r.norm() / h;
            if (q < 1) {
                std::pair<size_t, size_t> key(i, np->id);
                double d = gamma * springs[key];
                if (r.norm() > springs[key] + d) {
                    springs[key] += dt * alpha * (r.norm() - );
                }
            }
        }
    }
}


void sprint_displacements(std::vector<Particle> &pl, std::map<std::pair<size_t, size_t>, double> &springs);


void simulation_step(std::vector<Particle> &pl, Grid3D & grid) {
    retrieve_neiboroughs(pl, grid);
    apply_gravity(pl);
    viscosity_impulses(pl);
    avection(pl);
    double_density_relaxation(pl);
    update_vs(pl);
}

// void sprint_displacements(std::vector<Particle> &pl, std::vector<std::pair<size_t, size_t>> &springs) {

// }