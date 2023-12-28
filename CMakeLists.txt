cmake_minimum_required(VERSION 3.27)
project(black_scholes_MIMD_opt C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_FLAGS_RELEASE "-O0")
add_executable(black_scholes_MIMD_opt main.c)
