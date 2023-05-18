/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2023 Alexander Oleynichenko
 */

/*
 * Interface to the MyMathLib library:
 * http://www.mymathlib.com/
 */

#ifndef LIBGRPP_MYMATHLIB_GAMMA_H
#define LIBGRPP_MYMATHLIB_GAMMA_H

double Dawsons_Integral( double x );

double Gamma_Function(double x);

double Incomplete_Gamma_Function(double x, double nu);

double Entire_Incomplete_Gamma_Function(double x, double nu);

#endif // LIBGRPP_MYMATHLIB_GAMMA_H
