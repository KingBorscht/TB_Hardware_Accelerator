cmake_minimum_required(VERSION 3.16)
project(tb_mutation_caller LANGUAGES CXX)

add_executable(tb_mutation_caller mutation_caller.cpp)
target_compile_features(tb_mutation_caller PRIVATE cxx_std_17)
target_compile_options(tb_mutation_caller PRIVATE -Wall -Wextra -Wpedantic)
