/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2022 Alexander Oleynichenko
 */

#include "print_matrix.h"

#include <math.h>
#include <stdio.h>

void print_matrix(char *path, int dim, double *matrix)
{
    double const PRINT_THRESH = 1e-10;

    FILE *out_file = fopen(path, "w");

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            if (fabs(matrix[i * dim + j]) > PRINT_THRESH) {
                fprintf(out_file, "%4d%4d%30.16E\n", i + 1, j + 1, matrix[i * dim + j]);
            }
        }
    }

    fclose(out_file);
}


void print_matrix_lower_triangle(char *path, int dim, double *matrix)
{
    double const PRINT_THRESH = 1e-14;

    FILE *out_file = fopen(path, "w");

    fprintf(out_file, "%6d\n", dim);

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j <= i; j++) {
            if (fabs(matrix[i * dim + j]) > PRINT_THRESH) {
                fprintf(out_file, "%6d%6d%30.16E\n", i + 1, j + 1, matrix[i * dim + j]);
            }
        }
    }

    fclose(out_file);
}
