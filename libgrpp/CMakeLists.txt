#
#  libgrpp - a library for the evaluation of integrals over
#            generalized relativistic pseudopotentials.
#
#  Copyright (C) 2021-2023 Alexander Oleynichenko
#

cmake_minimum_required(VERSION 3.19)
project(libgrpp C)

set(CMAKE_C_STANDARD 11)

add_subdirectory(gsl)

add_compile_options(-O3)
#add_compile_options(-pg)
add_library(libgrpp
        angular_integrals.c
        binomial.c
        boys.c
        diff_gaussian.c
        factorial.c
        fortran.c
        full_grpp_integrals.c
        gfun.c
        grpp.c
        grpp_integrals_gradient.c
        gslbessel.c
        lmatrix.c
        norm_gaussian.c
        nuclear_attraction.c
        nuclear_models.c
        outercore_integrals.c
        overlap.c
        overlap_gradient.c
        parameters.c
        potential.c
        radial_type1_integral.c
        radial_type2_integral.c
        screening.c
        shell.c
        spherical_harmonics.c
        spin_orbit_integrals.c
        type1_integrals.c
        type1_mcmurchie_davidson.c
        type2_integrals.c
        utils.c
        )

target_link_libraries(libgrpp gsl -lm)
