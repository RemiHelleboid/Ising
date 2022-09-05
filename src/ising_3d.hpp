/**
 * @file ising_3d.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-08-30
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "ising_2d.hpp"

class ising_3d : public ising_base {
 private:
    std::size_t m_size_x;
    std::size_t m_size_y;
    std::size_t m_size_z;
    double      m_x_anisotropic_factor = 1.0;
    double      m_y_anisotropic_factor = 1.0;
    double      m_z_anisotropic_factor = 1.0;

 public:
    ising_3d(std::size_t size_x, std::size_t size_y, std::size_t size_z, double temperature = 1.0);

    void set_x_anisotropic_factor(double x_anisotropic_factor) { m_x_anisotropic_factor = x_anisotropic_factor; }
    void set_y_anisotropic_factor(double y_anisotropic_factor) { m_y_anisotropic_factor = y_anisotropic_factor; }
    void set_z_anisotropic_factor(double z_anisotropic_factor) { m_z_anisotropic_factor = z_anisotropic_factor; }

    double                                    get_spin(std::size_t x, std::size_t y, std::size_t z) const;
    void                                      set_spin(std::size_t x, std::size_t y, std::size_t z, double value);
    std::array<std::array<std::size_t, 3>, 8> get_neighbors(std::size_t x, std::size_t y, std::size_t z) const;

    double compute_energy(std::size_t x, std::size_t y, std::size_t z) const;
    double compute_total_energy() const override;
    double compute_specific_heat() const override;
    double compute_susceptibility() const override;

    void metropolis_step();

    ising_result metropolis_simulation(std::size_t nb_steps, const double convergence_threshold);
    void         metropolis_simulation_with_export(std::size_t nb_steps, const std::string& filename);

    void export_to_file(const std::string& filename) const;
};