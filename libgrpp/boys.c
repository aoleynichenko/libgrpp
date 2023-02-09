/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2023 Alexander Oleynichenko
 */

#include "boys.h"

#include <math.h>
#include <string.h>

#include "gslbessel.h"


/**
 * Boys function
 */
double boys(int n, double x)
{
    if (x > 1e-14) {
        double f = 2.0 * pow(x, n + 0.5);
        double g = gamma_function(n + 0.5);
        double p = incomplete_gamma_function(n + 0.5, x);

        return g * p / f;
    }
    else {
        return 1.0 / (n * 2.0 + 1.0);
    }
}


void boys_values(double x, int nmax, double *b)
{
    memset(b, 0, (nmax + 1) * sizeof(double));

    for (int n = 0; n <= nmax; n++) {
        b[n] = boys(n, x);
    }
}
