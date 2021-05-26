cmake_minimum_required(VERSION 3.18)
project(conditions)

set(CMAKE_CXX_STANDARD 20)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

add_subdirectory(lib)

# General

add_executable(dimension4sboxes dimension4sboxes.cc)
add_executable(dimension3experiments dimension3experiments.cc)
add_executable(dimension4experiments dimension4experiments.cc)
add_executable(dimension5experiments dimension5experiments.cc)
add_executable(dimension6experiments dimension6experiments.cc)
add_executable(dimension4timing dimension4timing.cc)
add_executable(dimension5timing dimension5timing.cc)
add_executable(dimension6timing dimension6timing.cc)
add_executable(dimension7timing dimension7timing.cc)

target_link_libraries(dimension4sboxes conditions_algorithms ntl gmp Threads::Threads)
target_link_libraries(dimension3experiments conditions_algorithms ntl gmp Threads::Threads)
target_link_libraries(dimension4experiments conditions_algorithms ntl gmp Threads::Threads)
target_link_libraries(dimension5experiments conditions_algorithms ntl gmp Threads::Threads)
target_link_libraries(dimension6experiments conditions_algorithms ntl gmp Threads::Threads)
target_link_libraries(dimension4timing conditions_algorithms ntl gmp Threads::Threads)
target_link_libraries(dimension5timing conditions_algorithms ntl gmp Threads::Threads)
target_link_libraries(dimension6timing conditions_algorithms ntl gmp Threads::Threads)
target_link_libraries(dimension7timing conditions_algorithms ntl gmp Threads::Threads)

# Cipher-specific

add_executable(trees_des trees_des.cc)
add_executable(trees_gift_f_delta trees_gift_f_delta.cc)
add_executable(trees_noekeon trees_noekeon.cc)
add_executable(trees_present trees_present.cc)
add_executable(trees_rectangle trees_rectangle.cc)
add_executable(trees_rectangle_f_delta trees_rectangle_f_delta.cc)
add_executable(trees_rectangle_fixed_components trees_rectangle_fixed_components.cc)
add_executable(trees_noekeon_fixed_components trees_noekeon_fixed_components.cc)
add_executable(trees_rectangle_multiple_outputs trees_rectangle_multiple_outputs.cc)
add_executable(trees_serpent trees_serpent.cc)

target_link_libraries(trees_des conditions_algorithms ntl gmp Threads::Threads)
target_link_libraries(trees_gift_f_delta conditions_algorithms ntl gmp Threads::Threads)
target_link_libraries(trees_noekeon conditions_algorithms ntl gmp Threads::Threads)
target_link_libraries(trees_present conditions_algorithms ntl gmp Threads::Threads)
target_link_libraries(trees_rectangle conditions_algorithms ntl gmp Threads::Threads)
target_link_libraries(trees_rectangle_f_delta conditions_algorithms ntl gmp Threads::Threads)
target_link_libraries(trees_rectangle_fixed_components conditions_algorithms ntl gmp Threads::Threads)
target_link_libraries(trees_noekeon_fixed_components conditions_algorithms ntl gmp Threads::Threads)
target_link_libraries(trees_rectangle_multiple_outputs conditions_algorithms ntl gmp Threads::Threads)
target_link_libraries(trees_serpent conditions_algorithms ntl gmp Threads::Threads)

# Linking omp
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(trees_serpent OpenMP::OpenMP_CXX)
endif()


