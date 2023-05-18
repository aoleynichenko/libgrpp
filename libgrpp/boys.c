/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2023 Alexander Oleynichenko
 */

#include "boys.h"

#include <math.h>
#include <string.h>

#include "scaled_mod_sph_bessel.h"
#include "mymathlib.h"


/**
 * Boys function
 */
double boys(int n, double x)
{
    if (x < 1e-14) {
        return 1.0 / (n * 2.0 + 1.0);
    }
    else if (x <= 0.01) { // Taylor expansion
        return 1.0 / (2.0 * n + 1.0)
               - x / (2.0 * n + 3.0)
               + 0.5 * x * x / (2.0 * n + 5.0)
               - x * x * x / ((2.0 * n + 7.0) * 6.0)
               + x * x * x * x / ((2.0 * n + 9.0) * 24.0)
               - x * x * x * x * x / ((2.0 * n + 11.0) * 120.0)
               + x * x * x * x * x * x / ((2.0 * n + 13.0) * 720.0);
    }
    else {
        double f = 2.0 * pow(x, n + 0.5);
        double g = Gamma_Function(n + 0.5);
        double p = Entire_Incomplete_Gamma_Function(x, n + 0.5);

        // (old GSL version:)
        //double g = gamma_function(n + 0.5);
        //double p = incomplete_gamma_function(n + 0.5, x);

        return g * p / f;
    }
}


void boys_values(double x, int nmax, double *b)
{
    memset(b, 0, (nmax + 1) * sizeof(double));

    for (int n = 0; n <= nmax; n++) {
        b[n] = boys(n, x);
    }
}
