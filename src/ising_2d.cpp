/**
 * @file ising_2d.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-08-29
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "ising_2d.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>

/**
 * @brief Construct a new ising 2d::ising 2d object.
 *
 * @param size_x
 * @param size_y
 * @param temperature
 */
ising_2d::ising_2d(std::size_t size_x, std::size_t size_y, double temperature)
    : m_size_x(size_x),
      m_size_y(size_y),
      m_spin_values(size_x * size_y),
      m_random_engine(std::random_device{}()),
      m_temperature(temperature) {}

/**
 * @brief Get the neighbor of the spin at position (x, y).
 *
 * @param x
 * @param y
 * @return std::array<std::pair<std::size_t, std::size_t>, 4>
 */
std::array<std::pair<std::size_t, std::size_t>, 4> ising_2d::get_neighbors(std::size_t x, std::size_t y) const {
    std::array<std::pair<std::size_t, std::size_t>, 4> neighbors;
    neighbors[0] = {x, (y + 1) % m_size_y};
    neighbors[1] = {(x + 1) % m_size_x, y};
    neighbors[2] = {x, (y + m_size_y - 1) % m_size_y};
    neighbors[3] = {(x + m_size_x - 1) % m_size_x, y};
    return neighbors;
}

/**
 * @brief Get the spin value at position (x, y).
 *
 * @param x
 * @param y
 * @return double
 */
double ising_2d::get_spin(std::size_t x, std::size_t y) const { return m_spin_values[x + y * m_size_x]; }

/**
 * @brief Set the spin value at position (x, y).
 *
 * @param x
 * @param y
 * @param value
 */
void ising_2d::set_spin(std::size_t x, std::size_t y, double value) { m_spin_values[x + y * m_size_x] = value; }

/**
 * @brief Initialize the spins randomly.
 *
 * @param probability
 *
 */
void ising_2d::initialize_random(double probability) {
    std::uniform_real_distribution<> double_distribution(0.0, 1.0);
    std::generate(m_spin_values.begin(), m_spin_values.end(), [&]() {
        return double_distribution(m_random_engine) < probability ? 1.0 : -1.0;
    });
}

/**
 * @brief Compute the energy of the spin at position (x, y).
 *
 * @param x
 * @param y
 * @return double
 */
double ising_2d::compute_energy(std::size_t x, std::size_t y) const {
    auto   neighbors = get_neighbors(x, y);
    double energy    = 0.0;
    for (auto& neighbor : neighbors) {
        energy += -get_spin(neighbor.first, neighbor.second) * get_spin(x, y) * (neighbor.first != x ? m_x_anisotropic_factor : 1.0) *
                  (neighbor.second != y ? m_y_anisotropic_factor : 1.0);
    }
    return energy;
}

/**
 * @brief Compute the total energy of the system.
 *
 * @return double
 */
double ising_2d::compute_total_energy() const {
    double energy = 0.0;
    for (std::size_t x = 0; x < m_size_x; x++) {
        for (std::size_t y = 0; y < m_size_y; y++) {
            energy += compute_energy(x, y);
        }
    }
    return energy;
}

/**
 * @brief Compute the total magnetization of the system.
 *
 * @return double
 */
double ising_2d::compute_total_magnetization() const { return std::accumulate(m_spin_values.begin(), m_spin_values.end(), 0.0); }

/**
 * @brief Compute the total magnetization of the system.
 *
 * @return double
 */
double ising_2d::compute_specific_heat() const {
    double energy = compute_total_energy();
    return energy * energy / (m_size_x * m_size_y);
}

/**
 * @brief Compute the total magnetization of the system.
 *
 * @return double
 */
double ising_2d::compute_susceptibility() const {
    double magnetization = compute_total_magnetization();
    return magnetization * magnetization / (m_size_x * m_size_y);
}

void ising_2d::metropolis_step() {
    m_number_modified_spins = 0;
    std::uniform_real_distribution<>           double_distribution(0.0, 1.0);
    std::uniform_int_distribution<std::size_t> int_distribution(0, m_size_y - 1);
    for (std::size_t idx_x = 0; idx_x < m_size_x; idx_x++) {
        for (std::size_t idx_y = 0; idx_y < m_size_y; idx_y++) {
            std::size_t x            = int_distribution(m_random_engine);
            std::size_t y            = int_distribution(m_random_engine);
            double      delta_energy = -2 * compute_energy(x, y);
            if (delta_energy <= 0.0) {
                set_spin(x, y, -get_spin(x, y));
                m_number_modified_spins++;
            } else {
                if (double_distribution(m_random_engine) < std::exp(-delta_energy / m_temperature)) {
                    set_spin(x, y, -get_spin(x, y));
                    m_number_modified_spins++;
                }
            }
        }
    }
}

void ising_2d::metropolis_step_parallel(int num_treads) {
    std::uniform_real_distribution<>           double_distribution(0.0, 1.0);
    std::uniform_int_distribution<std::size_t> int_distribution(0, m_size_y - 1);
#pragma omp parallel for
    for (std::size_t idx_x = 0; idx_x < m_size_x; idx_x++) {
        for (std::size_t idx_y = 0; idx_y < m_size_y; idx_y++) {
            std::size_t x            = int_distribution(m_random_engine);
            std::size_t y            = int_distribution(m_random_engine);
            double      delta_energy = -2 * compute_energy(x, y);
#pragma omp critical
            if (delta_energy <= 0.0) {
                set_spin(x, y, -get_spin(x, y));
            } else {
                if (double_distribution(m_random_engine) < std::exp(-delta_energy / m_temperature)) {
                    set_spin(x, y, -get_spin(x, y));
                }
            }
        }
    }
}

ising_result ising_2d::metropolis_simulation(std::size_t nb_steps, const double convergence_threshold) {
    double energy = compute_total_energy();

    for (std::size_t index_simulation = 0; index_simulation < nb_steps; index_simulation++) {
        metropolis_step();
        m_number_iterations++;
        double       ratio_modified_spins = static_cast<double>(m_number_modified_spins) / (m_size_x * m_size_y);
        double       new_energy           = compute_total_energy();
        ising_result result{compute_total_energy(), compute_total_magnetization(), compute_specific_heat(), compute_susceptibility()};
        if (ratio_modified_spins < convergence_threshold || (std::abs((new_energy - energy) / energy)) < convergence_threshold) {
            return result;
        }
        energy = new_energy;
    }
    ising_result result{compute_total_energy(), compute_total_magnetization(), compute_specific_heat(), compute_susceptibility()};
    return result;
}

void ising_2d::metropolis_simulation_with_export(std::size_t nb_steps, const std::string& filename) {
    std::ofstream file(filename + ".csv");
    file << "temperature,total_energy,total_magnetization,specific_heat,susceptibility" << std::endl;
    for (std::size_t index_simulation = 0; index_simulation < nb_steps; index_simulation++) {
        metropolis_step();

        std::ostringstream ss;
        ss << std::setw(5) << std::setfill('0') << index_simulation;
        std::string       str_index_simulation = ss.str();
        const std::string filename_iter        = filename + "_" + str_index_simulation + ".csv";
        export_to_file(filename_iter);
        file << m_temperature << "," << compute_total_energy() << "," << compute_total_magnetization() << "," << compute_specific_heat()
             << "," << compute_susceptibility() << std::endl;
        std::cout << "\r Iteration " << index_simulation << " / " << nb_steps << " (" << (index_simulation * 100.0 / nb_steps) << "%)"
                  << std::flush;
    }
}

void ising_2d::export_to_file(const std::string& filename) const {
    // std::cout << "Exporting to file " << filename << std::endl;
    std::ofstream file(filename);
    const double  size_x = static_cast<double>(m_size_x - 1);
    const double  size_y = static_cast<double>(m_size_y - 1);
    file << "X,Y,Spin\n";
    for (std::size_t x = 0; x < m_size_x; x++) {
        for (std::size_t y = 0; y < m_size_y; y++) {
            file << x / size_x << "," << y / size_y << "," << get_spin(x, y) << "\n";
        }
    }
    file.close();
}
