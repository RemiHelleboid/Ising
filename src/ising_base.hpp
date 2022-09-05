/**
 * @file ising_base.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-09-05
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

struct ising_result {
    double energy;
    double magnetization;
    double specific_heat;
    double susceptibility;
};

/**
 * @brief Base class for simple Ising model implementation.
 *
 */
class ising_base {
 protected:
    std::mt19937        m_random_engine;
    double              m_temperature;
    std::vector<double> m_spins;

    std::size_t m_number_iterations     = 0;
    std::size_t m_number_modified_spins = 0;

 public:
    ising_base(double temperature) : m_random_engine(std::random_device{}()), m_temperature(temperature){};
    ising_base(double temperature, std::size_t nb_spins)
        : m_random_engine(std::random_device{}()),
          m_temperature(temperature),
          m_spins(nb_spins){};
    virtual ~ising_base(){};

    void reset_spins() {
        std::fill(m_spins.begin(), m_spins.end(), 1.0);
        m_number_iterations     = 0;
        m_number_modified_spins = 0;
    }
    void resize_spins(std::size_t size) {
        m_spins.resize(size);
        reset_spins();
    }
    void        initialize_random(double probability);
    void        set_temperature(double temperature) { m_temperature = temperature; }
    std::size_t get_number_iterations() const { return m_number_iterations; }

    double         compute_total_magnetization() const;
    virtual double compute_total_energy() const   = 0;
    virtual double compute_specific_heat() const  = 0;
    virtual double compute_susceptibility() const = 0;
};
