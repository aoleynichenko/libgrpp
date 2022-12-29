/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2022 Alexander Oleynichenko
 */

#ifndef TEST_LIBGRPP_F90_X_OVERLAP_GRADIENTS_H
#define TEST_LIBGRPP_F90_X_OVERLAP_GRADIENTS_H

#include "shell_list.h"

void evaluate_overlap_integrals_gradients(
        int num_shells,
        libgrpp_shell_t **shell_list,
        molecule_t *molecule
);

#endif /* TEST_LIBGRPP_F90_X_OVERLAP_GRADIENTS_H */
