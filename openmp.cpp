#include "common.h"
#include <cmath>
#include <vector>
#include <omp.h>


// Apply the force from neighbor to particle
void apply_force(particle_t& particle, particle_t& neighbor) {
    // Calculate Distance
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

// Integrate the ODE
void move(particle_t& p, double size) {
    
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;
    // Bounce from walls
    while (p.x < 0 || p.x > size) {
        p.x = p.x < 0 ? -p.x : 2 * size - p.x;
        p.vx = -p.vx;
    }

    while (p.y < 0 || p.y > size) {
        p.y = p.y < 0 ? -p.y : 2 * size - p.y;
        p.vy = -p.vy;
    }
}

// 
void init_simulation(particle_t* parts, int num_parts, double size) {
    // No initialization needed for this optimized serial version.
}

void simulate_one_step(particle_t* parts, int num_parts, double size) {
    // Compute Forces - Optimized for O(n) interactions

    // 1. Divide the space into cells
    int num_cells_x = ceil(size / cutoff);
    int num_cells_y = ceil(size / cutoff);
    
    #pragma omp barrier

    std::vector<std::vector<std::vector<int>>> cells(num_cells_x, std::vector<std::vector<int>>(num_cells_y));
    // 2. Populate the cells
    //#pragma omp barrier// #pragma omp critical
    #pragma omp for nowait
    for (int i = 0; i < num_parts; ++i) {
        int cell_x = floor(parts[i].x / cutoff);
        int cell_y = floor(parts[i].y / cutoff);
        cells[cell_x][cell_y].push_back(i);
    }

    // 3. Iterate through cells and their neighbors
    //#pragma omp barrier
    // Each thread has its own local force arrays
    std::vector<double> local_ax(num_parts, 0.0);
    std::vector<double> local_ay(num_parts, 0.0);

    #pragma omp for nowait
    for (int i = 0; i < num_parts; ++i) {
        parts[i].ax = parts[i].ay = 0;
        int cell_x = parts[i].x / cutoff;
        int cell_y = parts[i].y / cutoff;

        for (int dx = -1; dx <= 1; ++dx) {
          
            for (int dy = -1; dy <= 1; ++dy) {
                int nx = cell_x + dx;
                int ny = cell_y + dy;

                if (nx >= 0 && nx < num_cells_x && ny >= 0 && ny < num_cells_y) {
		
                    for (int j : cells[nx][ny]) {
                        if (i != j) {
			    // Accumulate forces in thread-local storage
                                double dx = parts[j].x - parts[i].x;
                                double dy = parts[j].y - parts[i].y;
                                double r2 = dx * dx + dy * dy;

                                if (r2 > cutoff * cutoff) continue;

                                r2 = fmax(r2, min_r * min_r);
                                double r = sqrt(r2);
                                double coef = (1 - cutoff / r) / r2 / mass;

                                local_ax[i] += coef * dx;
                                local_ay[i] += coef * dy;
                        }
                    }
                }
            }
        }
    }

    // Move Particles
    // #pragma omp barrier
    #pragma omp for nowait
    for (int i = 0; i < num_parts; ++i) {
	parts[i].ax += local_ax[i];
        parts[i].ay += local_ay[i];
        move(parts[i], size);
    }
}
