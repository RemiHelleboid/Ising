file(GLOB ISING_SRC *.cpp)
file(GLOB ISING_INC *.hpp)

add_library(libising STATIC ${ISING_SRC} ${ISING_INC})
target_include_directories(libising PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} OpenMP::OpenMP_CXX)
