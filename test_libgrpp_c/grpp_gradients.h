/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2023 Alexander Oleynichenko
 */

#ifndef TEST_LIBGRPP_F90_X_GRPP_GRADIENTS_H
#define TEST_LIBGRPP_F90_X_GRPP_GRADIENTS_H

#include "abs_time.h"
#include "rpp.h"
#include "shell_list.h"

void evaluate_grpp_integrals_gradients(
        int num_shells,
        libgrpp_shell_t **shell_list,
        molecule_t *molecule,
        libgrpp_grpp_t **grpp_list
);

#endif /* TEST_LIBGRPP_F90_X_GRPP_GRADIENTS_H */
