/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2022 Alexander Oleynichenko
 */

/**
 * Calculation of nuclear attraction integrals.
 *
 * For the point charge nuclear model the recursive Obara-Saika scheme is used
 * to calculate nuclear attraction integrals. For details, see
 * T. Helgaker, P. Jorgensen, J. Olsen, Molecular Electronic-Structure Theory,
 * John Wiley & Sons Ltd, 2000.
 * Chapter 9.10.1, "The Obara-Saika scheme for one-electron Coulomb integrals"
 *
 * For the other three models,
 * - uniformly charged ball
 * - Gaussian nucleus
 * - Fermi nucleus,
 * the scheme is actually the same as for the type 1 (radially-local) ECP integrals.
 * Electrostatic potential V(r) induced by the finite nuclear charge distribution
 * is integrated numerically on a radial grid.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "boys.h"
#include "norm_gaussian.h"
#include "nuclear_models.h"
#include "libgrpp.h"
#include "utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct rpp_type1_data {
    double a;
    double b;
    double p;
    double mu;
    double A[3];
    double B[3];
    double C[3];
    double P[3];
    double R_QC_2;
    double K_abc;
};

double boys(int n, double x);

double evaluate_nuclear_attraction_integral_point_charge_primitive(
        double *charge_origin, int Z,
        double *origin_A, int n_A, int l_A, int m_A, double alpha_A,
        double *origin_B, int n_B, int l_B, int m_B, double alpha_B);

double nucattr_theta(struct rpp_type1_data *data, int N, int *ijklmn);

void evaluate_radially_local_potential_integral_primitive_gaussians(
        double *A, int n_cart_A, int *cart_list_A, double alpha_A,
        double *B, int n_cart_B, int *cart_list_B, double alpha_B,
        double *C, double (*potential)(double r, void *params),
        void *potential_params,
        double *matrix
);

double wrapper_coulomb_potential_point(double r, void *params);

double wrapper_coulomb_potential_ball(double r, void *params);

double wrapper_coulomb_potential_gaussian(double r, void *params);

double wrapper_coulomb_potential_fermi(double r, void *params);

double wrapper_coulomb_potential_fermi_bubble(double r, void *params);


/**
 * Calculates nuclear attraction integral between two shells
 * represented by contracted Gaussian functions.
 */
void libgrpp_nuclear_attraction_integrals(
        libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B,
        double *charge_origin, int charge, int nuclear_model, double *model_params,
        double *coulomb_matrix)
{
    int size_A = libgrpp_get_shell_size(shell_A);
    int size_B = libgrpp_get_shell_size(shell_B);

    double *buf = calloc(size_A * size_B, sizeof(double));

    memset(coulomb_matrix, 0, size_A * size_B * sizeof(double));

    // loop over primitives in contractions
    for (int i = 0; i < shell_A->num_primitives; i++) {
        for (int j = 0; j < shell_B->num_primitives; j++) {
            double coef_A_i = shell_A->coeffs[i];
            double coef_B_j = shell_B->coeffs[j];

            if (nuclear_model == LIBGRPP_NUCLEAR_MODEL_POINT_CHARGE) {

                // loop over cartesian functions inside the shells
                for (int m = 0; m < size_A; m++) {
                    for (int n = 0; n < size_B; n++) {
                        int n_A = shell_A->cart_list[3 * m + 0];
                        int l_A = shell_A->cart_list[3 * m + 1];
                        int m_A = shell_A->cart_list[3 * m + 2];
                        int n_B = shell_B->cart_list[3 * n + 0];
                        int l_B = shell_B->cart_list[3 * n + 1];
                        int m_B = shell_B->cart_list[3 * n + 2];

                        double s = evaluate_nuclear_attraction_integral_point_charge_primitive(
                                charge_origin, charge,
                                shell_A->origin, n_A, l_A, m_A, shell_A->alpha[i],
                                shell_B->origin, n_B, l_B, m_B, shell_B->alpha[j]
                        );

                        buf[m * size_B + n] = s;
                    }
                }
            }
            else if (nuclear_model == LIBGRPP_NUCLEAR_MODEL_CHARGED_BALL ||
                     nuclear_model == LIBGRPP_NUCLEAR_MODEL_GAUSSIAN ||
                     nuclear_model == LIBGRPP_NUCLEAR_MODEL_FERMI ||
                     nuclear_model == LIBGRPP_NUCLEAR_MODEL_FERMI_BUBBLE) {

                double params[10];
                params[0] = charge;

                /*
                 * choose nuclear model
                 */
                double (*charge_distrib_function)(double, void *) = NULL;
                if (nuclear_model == LIBGRPP_NUCLEAR_MODEL_CHARGED_BALL) {
                    params[1] = model_params[0]; // R_rms
                    charge_distrib_function = wrapper_coulomb_potential_ball;
                }
                else if (nuclear_model == LIBGRPP_NUCLEAR_MODEL_GAUSSIAN) {
                    params[1] = model_params[0]; // R_rms
                    charge_distrib_function = wrapper_coulomb_potential_gaussian;
                }
                else if (nuclear_model == LIBGRPP_NUCLEAR_MODEL_FERMI) {
                    params[1] = model_params[0]; // c
                    params[2] = model_params[1]; // a
                    charge_distrib_function = wrapper_coulomb_potential_fermi;
                }
                else {
                    params[1] = model_params[0]; // c
                    params[2] = model_params[1]; // a
                    params[3] = model_params[2]; // k
                    charge_distrib_function = wrapper_coulomb_potential_fermi_bubble;
                }

                /*
                 * calculate integrals for the shell pair
                 */
                evaluate_radially_local_potential_integral_primitive_gaussians(
                        shell_A->origin, size_A, shell_A->cart_list, shell_A->alpha[i],
                        shell_B->origin, size_B, shell_B->cart_list, shell_B->alpha[j],
                        charge_origin, charge_distrib_function, params, buf
                );
            }
            else {
                printf("LIBGRPP: unknown finite nuclear charge distribution model!\n");
                exit(0);
            }

            libgrpp_daxpy(size_A * size_B, coef_A_i * coef_B_j, buf, coulomb_matrix);
        }
    }

    free(buf);
}


/**
 * Three-dimensional overlap integral between two primitive Gaussian functions
 * g(r) = x^n y^l z^m exp(-alpha * r^2)
 */
double evaluate_nuclear_attraction_integral_point_charge_primitive(
        double *charge_origin, int Z,
        double *origin_A, int n_A, int l_A, int m_A, double alpha_A,
        double *origin_B, int n_B, int l_B, int m_B, double alpha_B)
{
    double N_A = gaussian_norm_factor(n_A, l_A, m_A, alpha_A);
    double N_B = gaussian_norm_factor(n_B, l_B, m_B, alpha_B);

    /*
     * auxiliary constants
     */
    struct rpp_type1_data data;
    data.a = alpha_A;
    data.b = alpha_B;
    data.p = alpha_A + alpha_B;
    data.mu = alpha_A * alpha_B / (alpha_A + alpha_B);
    for (int i = 0; i < 3; i++) {
        data.A[i] = origin_A[i];
        data.B[i] = origin_B[i];
        data.C[i] = charge_origin[i];
        data.P[i] = (alpha_A * origin_A[i] + alpha_B * origin_B[i]) / (alpha_A + alpha_B);
    }
    data.R_QC_2 = distance_squared(data.P, data.C);
    data.K_abc = exp(-data.mu * distance_squared(data.A, data.B));

    /*
     * start recursion
     */
    int ijklmn[6];
    ijklmn[0] = n_A;
    ijklmn[1] = n_B;
    ijklmn[2] = l_A;
    ijklmn[3] = l_B;
    ijklmn[4] = m_A;
    ijklmn[5] = m_B;

    return (-1) * Z * N_A * N_B * nucattr_theta(&data, 0, ijklmn);
}


/*
 * Obara-Saika recurrence relation for the Theta(N,i,j,k,l,m,n) function
 * (see Helgaker et al, Chapter 9.10.1)
 */
double nucattr_theta(struct rpp_type1_data *data, int N, int *ijklmn)
{
    /*
     * the base case of recursion
     */
    if (ijklmn[0] == 0 && ijklmn[1] == 0 &&
        ijklmn[2] == 0 && ijklmn[3] == 0 &&
        ijklmn[4] == 0 && ijklmn[5] == 0) {
        return 2 * M_PI / data->p * data->K_abc * boys(N, data->p * data->R_QC_2);
    }

    /*
     * recursion is performed first for the (i,j) indices (X),
     * then for (k,l) (Y), then for (l,m) (Z)
     */
    int coord;
    int i, j;
    if (ijklmn[0] != 0 || ijklmn[1] != 0) {
        coord = 0;
        i = ijklmn[0];
        j = ijklmn[1];
    }
    else if (ijklmn[2] != 0 || ijklmn[3] != 0) {
        coord = 1;
        i = ijklmn[2];
        j = ijklmn[3];
    }
    else {
        coord = 2;
        i = ijklmn[4];
        j = ijklmn[5];
    }

    /*
     * X_PA, X_PB, X_PC distances
     */
    double X_PA = data->P[coord] - data->A[coord];
    double X_PB = data->P[coord] - data->B[coord];
    double X_PC = data->P[coord] - data->C[coord];

    /*
     * downward recursion
     */
    double result = 0.0;
    double p = data->p;

    if (i >= j) { // downward step by i
        i -= 1;
        ijklmn[2 * coord] -= 1;
        result += X_PA * nucattr_theta(data, N, ijklmn);
        result -= X_PC * nucattr_theta(data, N + 1, ijklmn);
        if (i != 0) {
            ijklmn[2 * coord] -= 1;
            result += 0.5 / p * i * nucattr_theta(data, N, ijklmn);
            result -= 0.5 / p * i * nucattr_theta(data, N + 1, ijklmn);
            ijklmn[2 * coord] += 1;
        }
        if (j != 0) {
            ijklmn[2 * coord + 1] -= 1;
            result += 0.5 / p * j * nucattr_theta(data, N, ijklmn);
            result -= 0.5 / p * j * nucattr_theta(data, N + 1, ijklmn);
            ijklmn[2 * coord + 1] += 1;
        }
        ijklmn[2 * coord] += 1;
    }
    else { // downward step by j
        j -= 1;
        ijklmn[2 * coord + 1] -= 1;
        result += X_PB * nucattr_theta(data, N, ijklmn);
        result -= X_PC * nucattr_theta(data, N + 1, ijklmn);
        if (i != 0) {
            ijklmn[2 * coord] -= 1;
            result += 0.5 / p * i * nucattr_theta(data, N, ijklmn);
            result -= 0.5 / p * i * nucattr_theta(data, N + 1, ijklmn);
            ijklmn[2 * coord] += 1;
        }
        if (j != 0) {
            ijklmn[2 * coord + 1] -= 1;
            result += 0.5 / p * j * nucattr_theta(data, N, ijklmn);
            result -= 0.5 / p * j * nucattr_theta(data, N + 1, ijklmn);
            ijklmn[2 * coord + 1] += 1;
        }
        ijklmn[2 * coord + 1] += 1;
    }

    return result;
}


/**
 * wrappers for charge distribution functions.
 * are used to provide a unified interface to radially-local potentials.
 * the 'params' argument is unpacked, then the specific routines are invoked.
 */

double wrapper_coulomb_potential_ball(double r, void *params)
{
    double Z = ((double *) params)[0];
    double R_rms = ((double *) params)[1];

    return libgrpp_coulomb_potential_ball(r, Z, R_rms);
}


double wrapper_coulomb_potential_gaussian(double r, void *params)
{
    double Z = ((double *) params)[0];
    double R_rms = ((double *) params)[1];

    return libgrpp_coulomb_potential_gaussian(r, Z, R_rms);
}


double wrapper_coulomb_potential_fermi(double r, void *params)
{
    double Z = ((double *) params)[0];
    double c = ((double *) params)[1];
    double a = ((double *) params)[2];

    return libgrpp_coulomb_potential_fermi(r, Z, c, a);
}


double wrapper_coulomb_potential_fermi_bubble(double r, void *params)
{
    double Z = ((double *) params)[0];
    double c = ((double *) params)[1];
    double a = ((double *) params)[2];
    double k = ((double *) params)[3];

    return libgrpp_coulomb_potential_fermi_bubble(r, Z, c, a, k);
}
