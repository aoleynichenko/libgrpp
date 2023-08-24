/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2023 Alexander Oleynichenko
 */

#ifndef LIBGRPP_MOMENTUM_H
#define LIBGRPP_MOMENTUM_H

#include "libgrpp.h"

void libgrpp_momentum_integrals(libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B,
                                double *momentum_x_matrix, double *momentum_y_matrix, double *momentum_z_matrix);

#endif // LIBGRPP_MOMENTUM_H
