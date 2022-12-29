/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2022 Alexander Oleynichenko
 */

#ifndef EVAL_INTEGRALS_H_INCLUDED
#define EVAL_INTEGRALS_H_INCLUDED

#include "../libgrpp/libgrpp.h"
#include "shell_list.h"
#include "rpp.h"
#include "basis.h"
#include "xyz.h"

void evaluate_grpp_integrals(int num_shells, libgrpp_shell_t **shell_list,
                             molecule_t *molecule, grpp_t **grpp_list,
                             double *arep_matrix, double *so_x_matrix, double *so_y_matrix, double *so_z_matrix);

void evaluate_nuclear_attraction_integrals(int num_shells, libgrpp_shell_t **shell_list, molecule_t *molecule,
                                           double *nucattr_matrix, int nuclear_model);


#endif /* EVAL_INTEGRALS_H_INCLUDED */
