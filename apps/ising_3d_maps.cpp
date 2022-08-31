/**
 * @file ising_3d_maps.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-08-31
 *
 * @copyright Copyright (c) 2022
 *
 */
#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "ising_2d.hpp"
#include "ising_3d.hpp"

int main(int argc, char* argv[]) {
    std::size_t size_x               = 150;
    std::size_t size_y               = 150;
    std::size_t size_z               = 150;
    std::size_t nb_steps             = 100000;
    double      temperature          = 0.1;
    std::string filename             = "ising_2d_map";
    std::string out_dir              = ".";
    double      x_anisotropic_factor = 1.0;
    double      y_anisotropic_factor = 1.0;
    double      z_anisotropic_factor = 1.0;

    std::cout << "Usage: " << argv[0] << " [size_x] [size_y] [size_z] [nb_steps] [temperature] [filename] [outdir]" << std::endl;
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
        nb_steps = std::stoi(argv[4]);
    }
    if (argc > 5) {
        temperature = std::stod(argv[5]);
    }
    if (argc > 6) {
        filename = argv[6];
    }
    if (argc > 7) {
        out_dir = argv[7];
    }
    if (argc > 8) {
        x_anisotropic_factor = std::stod(argv[8]);
    }
    if (argc > 9) {
        y_anisotropic_factor = std::stod(argv[9]);
    }
    if (argc > 10) {
        z_anisotropic_factor = std::stod(argv[10]);
    }

    std::filesystem::create_directories(out_dir);
    ising_3d my_ising_3d(size_x, size_y, size_y, temperature);
    my_ising_3d.set_x_anisotropic_factor(x_anisotropic_factor);
    my_ising_3d.set_y_anisotropic_factor(y_anisotropic_factor);
    my_ising_3d.set_z_anisotropic_factor(z_anisotropic_factor);
    my_ising_3d.initialize_random(0.45);
    my_ising_3d.metropolis_simulation_with_export(nb_steps, out_dir + "/" + filename);

    return 0;
}