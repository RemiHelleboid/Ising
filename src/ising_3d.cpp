/**
 * @file ising_3d.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-08-30
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "ising_3d.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>

#include "ising_3d.hpp"

/**
 * @brief Construct a new ising 2d::ising 2d object.
 *
 * @param size_x
 * @param size_y
 * @param size_z
 * @param temperature
 */
ising_3d::ising_3d(std::size_t size_x, std::size_t size_y, std::size_t size_z, double temperature)
    : m_size_x(size_x),
      m_size_y(size_y),
      m_size_z(size_z),
      m_spin_values(size_x * size_y * size_z),
      m_random_engine(std::random_device{}()),
      m_temperature(temperature) {}

/**
 * @brief Get the neighbor of the spin at position (x, y, z).
 *
 * @param x
 * @param y
 * @param z
 * @return std::array<std::pair<std::size_t, std::size_t>, 8>
 */
std::array<std::array<std::size_t, 3>, 8> ising_3d::get_neighbors(std::size_t x, std::size_t y, std::size_t z) const {
    std::array<std::array<std::size_t, 3>, 8> neighbors;
    neighbors[0] = {x, (y + 1) % m_size_y, z};
    neighbors[1] = {x, (y - 1 + m_size_y) % m_size_y, z};
    neighbors[2] = {(x + 1) % m_size_x, y, z};
    neighbors[3] = {(x - 1 + m_size_x) % m_size_x, y, z};
    neighbors[4] = {x, y, (z + 1) % m_size_z};
    neighbors[5] = {x, y, (z - 1 + m_size_z) % m_size_z};
    neighbors[6] = {(x + 1) % m_size_x, (y + 1) % m_size_y, z};
    neighbors[7] = {(x - 1 + m_size_x) % m_size_x, (y - 1 + m_size_y) % m_size_y, z};
    return neighbors;
}

/**
 * @brief Get the spin value at position (x, y).
 *
 * @param x
 * @param y
 * @return double
 */
double ising_3d::get_spin(std::size_t x, std::size_t y, std::size_t z) const {
    return m_spin_values[x + y * m_size_x + z * m_size_x * m_size_y];
}

/**
 * @brief Set the spin value at position (x, y).
 *
 * @param x
 * @param y
 * @param value
 */
void ising_3d::set_spin(std::size_t x, std::size_t y, std::size_t z, double value) {
    m_spin_values[x + y * m_size_x + z * m_size_x * m_size_y] = value;
}

/**
 * @brief Initialize the spins randomly.
 *
 * @param probability
 *
 */
void ising_3d::initialize_random(double probability) {
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
double ising_3d::compute_energy(std::size_t x, std::size_t y, std::size_t z) const {
    auto   neighbors = get_neighbors(x, y, z);
    double energy    = 0.0;
    for (auto& neighbor : neighbors) {
        energy += -get_spin(neighbor[0], neighbor[1], neighbor[2]) * get_spin(x, y, z) * (neighbor[0] != x ? m_x_anisotropic_factor : 1.0) *
                  (neighbor[1] != y ? m_y_anisotropic_factor : 1.0) * (neighbor[2] != z ? m_z_anisotropic_factor : 1.0);
    }
    return energy;
}

/**
 * @brief Compute the total energy of the system.
 *
 * @return double
 */
double ising_3d::compute_total_energy() const {
    double energy = 0.0;
    for (std::size_t x = 0; x < m_size_x; x++) {
        for (std::size_t y = 0; y < m_size_y; y++) {
            for (std::size_t z = 0; z < m_size_z; z++) {
                energy += compute_energy(x, y, z);
            }
        }
    }
    return energy;
}

/**
 * @brief Compute the total magnetization of the system.
 *
 * @return double
 */
double ising_3d::compute_total_magnetization() const { return std::accumulate(m_spin_values.begin(), m_spin_values.end(), 0.0); }

/**
 * @brief Compute the total magnetization of the system.
 *
 * @return double
 */
double ising_3d::compute_specific_heat() const {
    double energy = compute_total_energy();
    return energy * energy / (m_size_x * m_size_y * m_size_z);
}

/**
 * @brief Compute the total magnetization of the system.
 *
 * @return double
 */
double ising_3d::compute_susceptibility() const {
    double magnetization = compute_total_magnetization();
    return magnetization * magnetization / (m_size_x * m_size_y * m_size_z);
}

void ising_3d::metropolis_step() {
    std::uniform_real_distribution<>           double_distribution(0.0, 1.0);
    std::uniform_int_distribution<std::size_t> int_distribution_x(0, m_size_x - 1);
    std::uniform_int_distribution<std::size_t> int_distribution_y(0, m_size_y - 1);
    std::uniform_int_distribution<std::size_t> int_distribution_z(0, m_size_z - 1);
    for (std::size_t idx_x = 0; idx_x < m_size_x; idx_x++) {
        for (std::size_t idx_y = 0; idx_y < m_size_y; idx_y++) {
            std::size_t x            = int_distribution_x(m_random_engine);
            std::size_t y            = int_distribution_y(m_random_engine);
            std::size_t z            = int_distribution_z(m_random_engine);
            double      delta_energy = -2 * compute_energy(x, y, z);
            if (delta_energy <= 0.0) {
                set_spin(x, y, z, -get_spin(x, y, z));
            } else {
                if (double_distribution(m_random_engine) < std::exp(-delta_energy / m_temperature)) {
                    set_spin(x, y, z, -get_spin(x, y, z));
                }
            }
        }
    }
}

ising_result ising_3d::metropolis_simulation(std::size_t nb_steps, const double convergence_threshold) {
    double energy = compute_total_energy();
    for (std::size_t index_simulation = 0; index_simulation < nb_steps; index_simulation++) {
        metropolis_step();
    }
    ising_result result{compute_total_energy(), compute_total_magnetization(), compute_specific_heat(), compute_susceptibility()};
    return result;
}

void ising_3d::metropolis_simulation_with_export(std::size_t nb_steps, const std::string& filename) {
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

void ising_3d::export_to_file(const std::string& filename) const {
    // std::cout << "Exporting to file " << filename << std::endl;
    std::ofstream file(filename);
    const double  size_x = static_cast<double>(m_size_x - 1);
    const double  size_y = static_cast<double>(m_size_y - 1);
    const double  size_z = static_cast<double>(m_size_z - 1);
    file << "X,Y,Z,Spin\n";
    for (std::size_t x = 0; x < m_size_x; x++) {
        for (std::size_t y = 0; y < m_size_y; y++) {
            for (std::size_t z = 0; z < m_size_y; z++) {
                file << x / size_x << "," << y / size_y << "," << z / size_z << "," << get_spin(x, y, z) << "\n";
            }
        }
    }
    file.close();
}
