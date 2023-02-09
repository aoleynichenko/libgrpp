/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2023 Alexander Oleynichenko
 */

#ifndef PRINT_MATRIX_H_INCLUDED
#define PRINT_MATRIX_H_INCLUDED

void print_matrix(char *path, int dim, double *matrix);
void print_matrix_lower_triangle(char *path, int dim, double *matrix);

#endif /* PRINT_MATRIX_H_INCLUDED */
