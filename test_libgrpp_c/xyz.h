/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2022 Alexander Oleynichenko
 */

#ifndef XYZ_H_INCLUDED
#define XYZ_H_INCLUDED

#include <stdio.h>

#include "elements.h"

typedef struct {
    int n_atoms;
    int *charges;
    double *coord_x;
    double *coord_y;
    double *coord_z;
} molecule_t;

molecule_t *new_molecule(int n_atoms, int *charges, double *x, double *y, double *z);

void delete_molecule(molecule_t *molecule);

molecule_t *read_molecule(char *path);

void print_molecule(FILE *out_file, molecule_t *molecule);

#endif /* XYZ_H_INCLUDED */
