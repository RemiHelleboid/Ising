/**
 * @file ising_base.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-09-05
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "ising_base.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

void ising_base::initialize_random(double probability) {
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    std::generate(m_spins.begin(), m_spins.end(), [&]() { return distribution(m_random_engine) < probability ? 1.0 : -1.0; });
}

double ising_base::compute_total_magnetization() const {
    double total_magnetization = std::accumulate(m_spins.begin(), m_spins.end(), 0.0) / static_cast<double>(m_spins.size());
}