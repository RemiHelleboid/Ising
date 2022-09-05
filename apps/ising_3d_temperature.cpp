/**
 * @file ising_3d_temperature.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-09-05
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "ising_2d.hpp"
#include "ising_3d.hpp"

void ising_3d_span_temperature(std::size_t        size_x,
                               std::size_t        size_y,
                               std::size_t        size_z,
                               double             temperature_min,
                               double             temperature_max,
                               double             temperature_step,
                               const std::string& filename) {
    std::ofstream file(filename);
    file << std::setprecision(10);
    std::size_t         nb_temperatures = (temperature_max - temperature_min) / temperature_step + 1;
    std::vector<double> temperatures(nb_temperatures);
    std::vector<double> energies(nb_temperatures);
    std::vector<double> magnetizations(nb_temperatures);
    std::vector<double> specific_heat(nb_temperatures);
    std::vector<double> susceptibility(nb_temperatures);
    for (std::size_t i = 0; i < nb_temperatures; ++i) {
        temperatures[i] = temperature_min + i * temperature_step;
    }
    std::cout << "nb temperatures: " << nb_temperatures << std::endl;
    const double threshold         = 1e-8;
    std::size_t  temperature_count = 0;
#pragma omp parallel for schedule(dynamic) reduction(+ : temperature_count)
    for (std::size_t i = 0; i < nb_temperatures; ++i) {
        ising_3d ising(size_x, size_y, size_z, temperatures[i]);
        ising.initialize_random(0.1);
        const int    nb_iter          = 20'000;
        ising_result result           = ising.metropolis_simulation(nb_iter, threshold);
        energies[i]                   = result.energy;
        magnetizations[i]             = result.magnetization;
        specific_heat[i]              = result.specific_heat;
        susceptibility[i]             = result.susceptibility;
        std::size_t number_iterations = ising.get_number_iterations();
        temperature_count++;
#pragma omp critical
        std::cout << "temperature: " << std::fixed << std::setprecision(2) << temperatures[i] << " converged in " << number_iterations
                  << " iterations. (" << temperature_count << ") / " << nb_temperatures << std::endl;
    }
    file << "temperature, energy, magnetization, specific_heat, susceptibility" << std::endl;
    for (std::size_t i = 0; i < nb_temperatures; ++i) {
        file << temperatures[i] << "," << energies[i] << "," << magnetizations[i] << "," << specific_heat[i] << "," << susceptibility[i]
             << std::endl;
    }
}

int main(int argc, char* argv[]) {
    std::size_t size_x           = 150;
    std::size_t size_y           = 150;
    std::size_t size_z           = 150;
    double      min_temperature  = 0.1;
    double      max_temperature  = 1.0;
    double      temperature_step = 0.1;
    std::string filename         = "ising_2d.csv";
    std::cout << "Usage: " << argv[0] << " [size_x] [size_y] [size_z] [min_temperature] [max_temperature] [temperature_step] [filename]"
              << std::endl;
    if (argc > 1) {
        size_x = std::stoi(argv[1]);
    }
    if (argc > 2) {
        size_y = std::stoi(argv[2]);
    }
    if (argc > 3) {
        size_z = std::stoi(argv[3]);
    }
    if (argc > 4) {
        min_temperature = std::stod(argv[4]);
    }
    if (argc > 5) {
        max_temperature = std::stod(argv[5]);
    }
    if (argc > 6) {
        temperature_step = std::stod(argv[6]);
    }
    if (argc > 7) {
        filename = argv[7];
    } else {
        filename = "ising_3d_" + std::to_string(size_x) + "x" + std::to_string(size_y) + "z" + std::to_string(size_z) + ".csv";
    }
    ising_3d_span_temperature(size_x, size_y, size_z, min_temperature, max_temperature, temperature_step, filename);
    return 0;
}