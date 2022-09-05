add_executable(temperatureIsing2d ising_2d_temperature.cpp)
target_link_libraries(temperatureIsing2d PUBLIC libising OpenMP::OpenMP_CXX)

add_executable(mapIsing2d ising_2d_maps.cpp)
target_link_libraries(mapIsing2d PUBLIC libising OpenMP::OpenMP_CXX)

add_executable(mapIsing3d ising_3d_maps.cpp)
target_link_libraries(mapIsing3d PUBLIC libising OpenMP::OpenMP_CXX)

add_executable(temperatureIsing3d ising_3d_temperature.cpp)
target_link_libraries(temperatureIsing3d PUBLIC libising OpenMP::OpenMP_CXX)
