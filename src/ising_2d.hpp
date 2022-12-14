/**
 * @file ising_2d.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-08-29
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

#include "ising_base.hpp"

class ising_2d : public ising_base {
 private:
    std::size_t m_size_x;
    std::size_t m_size_y;

    double m_x_anisotropic_factor = 1.0;
    double m_y_anisotropic_factor = 1.0;

 public:
    ising_2d(std::size_t size_x, std::size_t size_y, double temperature = 1.0);

    void set_x_anisotropic_factor(double x_anisotropic_factor) { m_x_anisotropic_factor = x_anisotropic_factor; }
    void set_y_anisotropic_factor(double y_anisotropic_factor) { m_y_anisotropic_factor = y_anisotropic_factor; }

    double                                             get_spin(std::size_t x, std::size_t y) const;
    void                                               set_spin(std::size_t x, std::size_t y, double value);
    std::array<std::pair<std::size_t, std::size_t>, 4> get_neighbors(std::size_t x, std::size_t y) const;

    double compute_energy(std::size_t x, std::size_t y) const;
    double compute_total_energy() const;
    double compute_total_magnetization() const;
    double compute_specific_heat() const;
    double compute_susceptibility() const;

    void metropolis_step();
    void metropolis_step_parallel(int num_treads);

    ising_result metropolis_simulation(std::size_t nb_steps, const double convergence_threshold);
    void         metropolis_simulation_with_export(std::size_t nb_steps, const std::string& filename);

    void export_to_file(const std::string& filename) const;
};