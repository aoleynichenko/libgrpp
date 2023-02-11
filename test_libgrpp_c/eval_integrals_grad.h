/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2023 Alexander Oleynichenko
 */

#ifndef TEST_LIBGRPP_F90_X_EVAL_INTEGRALS_GRAD_H
#define TEST_LIBGRPP_F90_X_EVAL_INTEGRALS_GRAD_H

#include "../libgrpp/libgrpp.h"
#include "molecule.h"

void evaluate_grpp_integrals_gradient(
        int num_shells,
        libgrpp_shell_t **shell_list,
        molecule_t *molecule,
        libgrpp_grpp_t **grpp_list,
        double **gradient_arep,
        double **gradient_so_x,
        double **gradient_so_y,
        double **gradient_so_z
);

void evaluate_overlap_integrals_gradient(
        int num_shells,
        libgrpp_shell_t **shell_list,
        molecule_t *molecule,
        double **gradient
);

#endif //TEST_LIBGRPP_F90_X_EVAL_INTEGRALS_GRAD_H
