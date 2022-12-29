#
#  libgrpp - a library for the evaluation of integrals over
#            generalized relativistic pseudopotentials.
#
#  Copyright (C) 2021-2022 Alexander Oleynichenko
#

cmake_minimum_required(VERSION 3.19)
project(libgrpp C)

set(CMAKE_C_STANDARD 11)

add_subdirectory(gsl)

add_compile_options(-O3)
#add_compile_options(-pg)
add_library(libgrpp
        spherical_harmonics.c
        binomial.c
        boys.c
        gfun.c
        factorial.c
        norm_gaussian.c
        gslbessel.c
        type1_integrals.c
        type1_mcmurchie_davidson.c
        angular_integrals.c
        type2_integrals.c
        spin_orbit_integrals.c
        utils.c
        radial_type1_integral.c
        radial_type2_integral.c
        overlap.c
        nuclear_models.c
        nuclear_attraction.c
        outercore_integrals.c
        fortran.c
        parameters.c
        lmatrix.c
        potential.c
        shell.c
        screening.c
        )

target_link_libraries(libgrpp gsl -lm)
