/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2022 Alexander Oleynichenko
 */

#ifndef RPP_H_INCLUDED
#define RPP_H_INCLUDED

#include <stdio.h>

#include "../libgrpp/libgrpp.h"

#define MAX_NUM_OC_SHELLS 100

typedef struct {
    int n_arep;
    int n_esop;
    int n_oc_shells;
    libgrpp_potential_t *U_L;
    libgrpp_potential_t **U_arep;
    libgrpp_potential_t **U_esop;
    libgrpp_potential_t **U_oc;
    libgrpp_shell_t **oc_shells;
} grpp_t;

grpp_t *read_grpp(char *path, int nuc_charge);

void delete_grpp(grpp_t *grpp);

void print_grpp(FILE *out, grpp_t *grpp);

#endif /* RPP_H_INCLUDED */
