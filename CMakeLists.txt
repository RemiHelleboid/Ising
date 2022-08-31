cmake_minimum_required(VERSION 3.15...3.20)

project(
    Ising++
    VERSION "0.0.0"
    DESCRIPTION "Ising++ Project"
    LANGUAGES C CXX)

include(FetchContent)
include(CheckCXXCompilerFlag)
include(CMakePrintHelpers)

set(CMAKE_CXX_STANDARD 20)

add_compile_definitions(CMAKE_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

# add_compile_definitions(ARMIN_VERSION = $CMAKE_PROJECT_VERSION)

# this variable controls the default value of the options which are normally set
# to ON (useful if you want to configure a minimal version of Gmsh: e.g. "cmake
# -DDEFAULT=0 -DENABLE_POST=1 -DENABLE_PARSER=1")
set(DEFAULT
    ON
    CACHE INTERNAL "Default value for enabled-by-default options")

set(default_build_type "Release")

if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
    set(default_build_type "RelWithDebInfo")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(
        STATUS
        "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE
        "${default_build_type}"
        CACHE STRING "Choose the type of build." FORCE)

    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
        "MinSizeRel" "RelWithDebInfo")
endif()

if(ENABLE_CLANG_TIDY)
    message("Clang tidy at compilation is enabled")
    set(CMAKE_CXX_CLANG_TIDY
        clang-tidy -header-filter={CMAKE_SOURCE_DIR}/src/Mesh;
        -checks=*,readability*,modernize*,-modernize-use-trailing-return-type,-fuchsia-default-arguments-calls,-hicpp-use-emplace,-modernize-use-emplace)
endif()

find_package(OpenMP)

add_subdirectory(src)
add_subdirectory(apps)
