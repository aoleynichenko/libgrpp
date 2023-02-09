/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2023 Alexander Oleynichenko
 */

#ifndef SHELL_LIST_H_INCLUDED
#define SHELL_LIST_H_INCLUDED

#include "../libgrpp/libgrpp.h"
#include "molecule.h"
#include "basis.h"


libgrpp_shell_t **construct_shell_list(molecule_t *mol, basis_set_t **basis_sets, int *num_shells);

void delete_shell_list(libgrpp_shell_t **shell_list, int num_shells);

int calculate_basis_dim(libgrpp_shell_t **shell_list, int num_shells);

void print_shell_list(libgrpp_shell_t **shell_list, int num_shells);

#endif /* SHELL_LIST_H_INCLUDED */
