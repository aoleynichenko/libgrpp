/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2022 Alexander Oleynichenko
 */

#include "xyz.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


molecule_t *new_molecule(int n_atoms, int *charges, double *x, double *y, double *z)
{
    molecule_t *molecule = (molecule_t *) malloc(sizeof(molecule_t) * 1);

    molecule->charges = (int *) malloc(sizeof(int) * n_atoms);
    molecule->coord_x = (double *) malloc(sizeof(double) * n_atoms);
    molecule->coord_y = (double *) malloc(sizeof(double) * n_atoms);
    molecule->coord_z = (double *) malloc(sizeof(double) * n_atoms);

    molecule->n_atoms = n_atoms;
    for (int i = 0; i < n_atoms; i++) {
        molecule->charges[i] = charges[i];
        molecule->coord_x[i] = x[i];
        molecule->coord_y[i] = y[i];
        molecule->coord_z[i] = z[i];
    }

    return molecule;
}


void delete_molecule(molecule_t *molecule)
{
    free(molecule->charges);
    free(molecule->coord_x);
    free(molecule->coord_y);
    free(molecule->coord_z);
    free(molecule);
}


molecule_t *read_molecule(char *path)
{
    FILE *inp_file;
    int n_atoms;
    int *charges;
    double *x, *y, *z;
    char buf[256];
    molecule_t *molecule = NULL;

    // open xyz file
    inp_file = fopen(path, "r");
    if (inp_file == NULL) {
        printf("Cannot open xyz file '%s\n'", path);
        return NULL;
    }

    // read total number of atoms
    if (fscanf(inp_file, "%d", &n_atoms) != 1) {
        printf("Error while reading xyz file\n");
        return NULL;
    }

    // read xyz coordinates for each atom
    charges = (int *) malloc(sizeof(int) * n_atoms);
    x = (double *) malloc(sizeof(double) * n_atoms);
    y = (double *) malloc(sizeof(double) * n_atoms);
    z = (double *) malloc(sizeof(double) * n_atoms);

    for (int i = 0; i < n_atoms; i++) {
        if (fscanf(inp_file, "%s%lf%lf%lf", buf, &x[i], &y[i], &z[i]) != 4) {
            printf("Error while reading xyz file\n");
            molecule = NULL;
            goto cleanup;
        }

        int nuc_charge = get_element_nuc_charge(buf);
        if (nuc_charge == -1) {
            printf("Error while reading xyz file\n");
            printf("Unknown element: %s\n", buf);
            molecule = NULL;
            goto cleanup;
        }
        charges[i] = nuc_charge;
    }

    fclose(inp_file);

    molecule = new_molecule(n_atoms, charges, x, y, z);

    cleanup:
    free(charges);
    free(x);
    free(y);
    free(z);

    return molecule;
}


void print_molecule(FILE *out_file, molecule_t *molecule)
{
    fprintf(out_file, "\n\t\tgeometry:\n");
    fprintf(out_file, "\t\t---------\n");
    for (int i = 0; i < molecule->n_atoms; i++) {
        fprintf(out_file, "  %3d%12.8f%12.8f%12.8f\n", molecule->charges[i],
               molecule->coord_x[i], molecule->coord_y[i], molecule->coord_z[i]);
    }
    fprintf(out_file, "\n");
}
