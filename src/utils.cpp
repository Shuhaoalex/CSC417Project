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

void double_density_relaxation(std::vector<Particle> &pl, constants &c) {
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
        double P = c.k * (rho - c.rho0);
        double P_near = c.k_near * rho_near;
        Eigen::Vector3d dx = Eigen::Vector3d::Zero();
        for (auto& np: par.neighbors) {
            Eigen::Vector3d r = np->x - par.x;
            double q = r.norm() / h;
            if (q < 1) {
                Eigen::Vector3d D = pow(c.dt, 2) * (P * (1-q) + P_near * pow(1-q, 2)) / 2 * r.normalized();
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

void viscosity_impulses(std::vector<Particle> &pl, constants &c) {
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
                        Eigen::Vector3d I = c.dt * (1-q) * (c.sigma * u + c.beta * u * u) / 2 * r.normalized();
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

void apply_gravity(std::vector<Particle> &pl, constants &c) {
    #pragma omp parallel for
    for (int i = 0; i < pl.size(); ++i) {
        pl[i].v += c.dt * c.g;
    }
}

void avection(std::vector<Particle> &pl, constants &c) {
    #pragma omp parallel for
    for (int i = 0; i < pl.size(); ++i) {
        pl[i].x_prev = pl[i].x;
        pl[i].x += c.dt * pl[i].v;
    }
}

void update_vs(std::vector<Particle> &pl, constants &c) {
    #pragma omp parallel for
    for (int i = 0; i < pl.size(); ++i) {
        auto&p = pl[i];
        p.v = (p.x - p.x_prev) / c.dt;
        for (int ax = 0; ax < 3; ++ax) {
            if (p.x(ax) < 0) {
                p.x(ax) = 0;
                p.v(ax) *= -c.collision_e;
            } else if (p.x(ax) > box_size) {
                p.x(ax) = box_size;
                p.v(ax) *= -c.collision_e;
            }
        }
    }
}

void adjustSprings(std::vector<Particle> &pl, std::map<std::pair<size_t, size_t>, Spring> &springs, constants &c) {
    for (size_t i = 0; i < pl.size(); ++i) {
        auto& par = pl[i];
        for (auto& np: par.neighbors) {
            if (par.id < np->id) {
                double r = (np->x - par.x).norm();
                if (r < h) {
                    std::pair<size_t, size_t> key(i, np->id);
                    double Lij = springs[key].value;
                    double d = c.gamma * Lij;
                    if (r > Lij + d) {
                        springs[key].value += c.dt * c.alpha * (r - Lij - d);
                    } else if (r < Lij - d) {
                        springs[key].value -= c.dt * c.alpha * (Lij - d - r);
                    }
                }
            }
        }
    }
    auto it = springs.begin();
    while (it != springs.end()) {
        size_t i = it->first.first;
        size_t j = it->first.second;
        // if ((pl[j].x - pl[i].x).norm() >= h) {
        //     it = springs.erase(it);
        // }
        if (it->second.value >= h) {
            it = springs.erase(it);
        } else {
            ++it;
        }
    }
    // for (auto it = springs.begin(); it != springs.end(); ++it) {
    //     // size_t i = it->first.first;
    //     // size_t j = it->first.second;
    //     // if ((pl[j].x - pl[i].x).norm() >= h) {
    //     //     springs.erase(it);
    //     // }
        
    // }
}


void spring_displacements(std::vector<Particle> &pl, std::map<std::pair<size_t, size_t>, Spring> &springs, constants &c) {
    for (auto it = springs.begin(); it != springs.end(); ++it) {
        size_t i = it->first.first;
        size_t j = it->first.second;
        Eigen::Vector3d r = pl[j].x - pl[i].x;
        Eigen::Vector3d D  = c.dt * c.dt * c.k_spring * (1 - it->second.value / h) * (it->second.value - r.norm()) * r.normalized();
        pl[i].x -= D / 2;
        pl[j].x += D / 2;
    }
}


void simulation_step(std::vector<Particle> &pl, Grid3D & grid, std::map<std::pair<size_t, size_t>, Spring> &springs, constants &c) {
    retrieve_neiboroughs(pl, grid);
    apply_gravity(pl, c);
    viscosity_impulses(pl, c);
    avection(pl, c);
    // adjustSprings(pl, springs, c);
    // spring_displacements(pl, springs, c);
    double_density_relaxation(pl, c);
    update_vs(pl, c);
}

// void sprint_displacements(std::vector<Particle> &pl, std::vector<std::pair<size_t, size_t>> &springs) {

// }