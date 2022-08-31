/**
 * @file ising_2d_maps.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-08-30
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

int main(int argc, char* argv[]) {
    std::size_t size_x               = 150;
    std::size_t size_y               = 150;
    std::size_t nb_steps             = 100000;
    double      temperature          = 0.1;
    std::string filename             = "ising_2d_map";
    std::string out_dir              = ".";
    double      x_anisotropic_factor = 1.0;
    double      y_anisotropic_factor = 1.0;

    std::cout << "Usage: " << argv[0] << " [size_x] [size_y] [nb_steps] [temperature] [filename] [outdir]" << std::endl;
    if (argc > 1) {
        size_x = std::stoi(argv[1]);
    }
    if (argc > 2) {
        size_y = std::stoi(argv[2]);
    }
    if (argc > 3) {
        nb_steps = std::stoi(argv[3]);
    }
    if (argc > 4) {
        temperature = std::stod(argv[4]);
    }
    if (argc > 5) {
        filename = argv[5];
    }
    if (argc > 6) {
        out_dir = argv[6];
    }
    if (argc > 7) {
        x_anisotropic_factor = std::stod(argv[7]);
    }
    if (argc > 8) {
        y_anisotropic_factor = std::stod(argv[8]);
    }
    if (argc <= 6) {
        out_dir = "ising2d_results_" + std::to_string(size_x) + "x" + std::to_string(size_y) + "_T" + std::to_string(temperature) + "/";
    }
    std::cout << "Results stored in: " << out_dir << " folder.\n";

    std::filesystem::create_directories(out_dir);
    ising_2d my_ising_2d(size_x, size_y, temperature);
    my_ising_2d.set_x_anisotropic_factor(x_anisotropic_factor);
    my_ising_2d.set_y_anisotropic_factor(y_anisotropic_factor);
    my_ising_2d.initialize_random(0.45);
    my_ising_2d.metropolis_simulation_with_export(nb_steps, out_dir + "/" + filename);

    const std::string python_script = CMAKE_SOURCE_DIR + std::string("/python/parse_Ising2d.py");
    const std::string python_call = "python3 " + python_script + " -d " + out_dir + " -s 0";
    int python_success = system(python_call.c_str());

    return 0;
}