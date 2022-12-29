/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2022 Alexander Oleynichenko
 */

#ifndef LIBGRPP_GSLBESSEL_H
#define LIBGRPP_GSLBESSEL_H

double modified_bessel_scaled(int n, double x);

double gamma_function(double x);

double incomplete_gamma_function(double a, double x);

#endif //LIBGRPP_GSLBESSEL_H
