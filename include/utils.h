#ifndef __UTILS_H__
#define __UTILS_H__
#include <Eigen/Core>
#include "Grid3D.h"
#include <vector>
#include <utility>
void retrieve_neiboroughs(std::vector<Particle> &pl, Grid3D & grid);
// void adjustSprings(std::vector<Particle> &pl, std::vector<std::pair<size_t, size_t>> &springs);
// void sprint_displacements(std::vector<Particle> &pl, std::vector<std::pair<size_t, size_t>> &springs);
void double_density_relaxation(std::vector<Particle> &pl);
void viscosity_impulses(std::vector<Particle> &pl);
void apply_gravity(std::vector<Particle> &pl);
void avection(std::vector<Particle> &pl);
void update_vs(std::vector<Particle> &pl);
void simulation_step(std::vector<Particle> &pl, Grid3D & grid);
#endif