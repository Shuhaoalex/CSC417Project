#ifndef __GRID_H__
#define __GRID_H__
#include <Eigen/Core>
#include <vector>
#include <cmath>
#include <algorithm>
#include "Particle.h"
#include "LinkedList.h"
#include <iostream>

class Grid3D {
    private:
    std::vector<LinkedListNode<Particle *> *> grid_list;
    std::vector<LinkedListNode<Particle *>> particle_list;
    int sz;
    double dg;

    inline int flatten_idx(Eigen::Vector3i grid_loc) {
        return grid_loc(0) * sz * sz + grid_loc(1) * sz + grid_loc(2);
    }

    public:
    inline Grid3D(int size, double query_radius):sz(size), dg(query_radius){
        grid_list.resize(pow(sz, 3));
        for (int i = 0; i < grid_list.size(); ++i) {
            grid_list[i] = NULL;
        }
    }

    inline void clear() {
        for (int i = 0; i < grid_list.size(); ++i) {
            grid_list[i] = NULL;
        }
        particle_list.clear();   
    }

    inline void update_particles_pos() {
        for (int i = 0; i < particle_list.size(); ++i) {
            Eigen::Vector3d relative_coord = particle_list[i].payload->x / dg;
            Eigen::Vector3i grid_loc = relative_coord.array().floor().cast<int>().min(sz - 1);
            particle_list[i].payload->grid_loc = grid_loc;
            int list_idx = flatten_idx(grid_loc);
            particle_list[i].remove();
            particle_list[i].insert(grid_list[list_idx]);
        }
    }

    inline void insert_particles(std::vector<Particle> &pl) {
        for (int i = 0; i < pl.size(); ++i) {
            particle_list.emplace_back(&pl[i]);
        }
        update_particles_pos();
    }

    inline void get_neighbors(const Particle& p, std::vector<Particle*> &result) {
        Eigen::Vector3i p_loc = p.grid_loc;
        result.clear();
        for (int i = std::max(0, p_loc(0) - 1); i < std::min(sz, p_loc(0) + 2); ++i) {
            for (int j = std::max(0, p_loc(1) - 1); j < std::min(sz, p_loc(1) + 2); ++j) {
                for (int k = std::max(0, p_loc(2) - 1); k < std::min(sz, p_loc(2) + 2); ++k) {
                    Eigen::Vector3i curr_loc(i, j, k);
                    int list_idx = flatten_idx(curr_loc);
                    LinkedListNode<Particle *> * curr = grid_list[list_idx];
                    while (curr) {
                        if (p.id != curr->payload->id) {
                            if ((p.x - curr->payload->x).norm() < dg) {
                                result.emplace_back(curr->payload);
                            }
                        }
                        curr = curr->next();
                    }
                }
            }
        }
    }
};

#endif