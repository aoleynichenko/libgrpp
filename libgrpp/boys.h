/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2022 Alexander Oleynichenko
 */

#ifndef LIBGRPP_BOYS_H
#define LIBGRPP_BOYS_H

double boys(int n, double x);

void boys_values(double x, int nmax, double *b);

#endif // LIBGRPP_BOYS_H
