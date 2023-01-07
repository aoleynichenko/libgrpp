#
#  libgrpp - a library for the evaluation of integrals over
#            generalized relativistic pseudopotentials.
#
#  Copyright (C) 2021-2022 Alexander Oleynichenko
#

cmake_minimum_required(VERSION 3.19)
project(test_libgrpp_c.x C)
project(test_libgrpp_f90.x Fortran)

set(CMAKE_C_STANDARD 11)

add_subdirectory(libgrpp)

#add_compile_options(-pg)
add_compile_options(-O3)
add_executable(test_libgrpp_c.x
        test_libgrpp_c/main.c
        test_libgrpp_c/xyz.c
        test_libgrpp_c/basis.c
        test_libgrpp_c/rpp.c
        test_libgrpp_c/elements.c
        test_libgrpp_c/shell_list.c
        test_libgrpp_c/eval_integrals.c
        test_libgrpp_c/overlap_gradients.c
        test_libgrpp_c/grpp_gradients.c
        test_libgrpp_c/print_matrix.c
        test_libgrpp_c/abs_time.c
        )

add_executable(test_libgrpp_f90.x
        test_libgrpp_f90/main.f90
        test_libgrpp_f90/xyz.f90
        test_libgrpp_f90/basis.f90
        test_libgrpp_f90/rpp.f90
        test_libgrpp_f90/evalints.f90
        test_libgrpp_f90/libgrpp.f90
        )


target_link_libraries(test_libgrpp_c.x libgrpp -lm)
target_link_libraries(test_libgrpp_f90.x libgrpp -lm)


# testing binary
add_executable(libgrpp_compare_matrices.x test/compare_matrices.c)

# enable testing functionality
enable_testing()

# define tests

add_test(
        NAME Cs
        WORKING_DIRECTORY ../test/Cs
        COMMAND bash run_test.sh
)

add_test(
        NAME H2O
        WORKING_DIRECTORY ../test/H2O
        COMMAND bash run_test.sh
)

add_test(
        NAME ThO
        WORKING_DIRECTORY ../test/ThO
        COMMAND bash run_test.sh
)

add_test(
        NAME UO2
        WORKING_DIRECTORY ../test/UO2
        COMMAND bash run_test.sh
)


