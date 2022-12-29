/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2022 Alexander Oleynichenko
 */

#ifndef BASIS_H_INCLUDED
#define BASIS_H_INCLUDED

#include <stdio.h>

typedef struct {
    int L;
    int n_primitives;
    double *coeffs;
    double *alpha;
} basis_fun_t;

basis_fun_t *new_basis_function(int L, int n_primitives, double *coeffs, double *alpha);

void delete_basis_function(basis_fun_t *fun);


typedef struct {
    int n_fun;
    int capacity;
    basis_fun_t **fun_list;
} basis_set_t;

basis_set_t *new_basis_set();

void delete_basis_set(basis_set_t *set);

void add_basis_function(basis_set_t *set, basis_fun_t *fun);

basis_set_t *read_basis_set(char *path, int nuc_charge);

void print_basis_set(FILE *out, basis_set_t *set);

#endif /* BASIS_H_INCLUDED */
