/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2023 Alexander Oleynichenko
 */

/*
 * This example shown how to calculate numerical and analytic gradients
 * of overlap integrals using the LIBGRPP library.
 *
 * For details, see
 * T. Helgaker, P. Jorgensen, J. Olsen, Molecular Electronic-Structure Theory,
 * John Wiley & Sons Ltd, 2000.
 * Chapter 9.3.1, "Overlap integrals"
 */

#include "overlap_gradients.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "abs_time.h"
#include "shell_list.h"

void libgrpp_overlap_integrals_num_grad(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        double *nuc_coord,
        double **grad
);

void libgrpp_overlap_integrals_analytic_grad(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        double *nuc_coord,
        double **grad
);

void overlap_gradient_diff_bra_contribution(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        double **grad,
        double factor
);

void overlap_gradient_diff_bra_overlap_integrals(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        double **overlap_down,
        double **overlap_up,
        int *cart_size_down,
        int *cart_size_up
);

void differentiate_shell(libgrpp_shell_t *shell, libgrpp_shell_t **shell_minus, libgrpp_shell_t **shell_plus);

int nlm_to_linear(int *nlm);

int points_are_equal(double *a, double *b);

double norm_factor(double alpha, int L);

void print_gradients(libgrpp_shell_t *bra, libgrpp_shell_t *ket, double **grad);

double **alloc_gradients(libgrpp_shell_t *bra, libgrpp_shell_t *ket);

void dealloc_gradients(double **grad);


/**
 * This function evaluates gradients of overlap integrals with respect to coordinates
 * of all atoms in a given molecule.
 */
void evaluate_overlap_integrals_gradients(
        int num_shells,
        libgrpp_shell_t **shell_list,
        molecule_t *molecule
)
{
    double time_num = 0.0;
    double time_ana = 0.0;

    int ioffset = 0;
    for (int ishell = 0; ishell < num_shells; ishell++) {

        libgrpp_shell_t *bra = shell_list[ishell];
        int bra_dim = libgrpp_get_shell_size(bra);

        int joffset = 0;
        for (int jshell = 0; jshell < num_shells; jshell++) {

            libgrpp_shell_t *ket = shell_list[jshell];
            int ket_dim = libgrpp_get_shell_size(ket);

            if (points_are_equal(bra->origin, ket->origin)) {
                continue;
            }

            printf("\n\n");
            printf("Shell %4d   L = %1d    %10.6f%10.6f%10.6f\n", ishell + 1, bra->L, bra->origin[0], bra->origin[1],
                   bra->origin[2]);
            printf("Shell %4d   L = %1d    %10.6f%10.6f%10.6f\n", jshell + 1, ket->L, ket->origin[0], ket->origin[1],
                   ket->origin[2]);

            for (int iatom = 0; iatom < molecule->n_atoms; iatom++) {
                printf("\n");
                printf("wrt %6.3f%6.3f%6.3f:\n",
                       molecule->coord_x[iatom], molecule->coord_y[iatom], molecule->coord_z[iatom]);

                double diff_origin[3];
                diff_origin[0] = molecule->coord_x[iatom];
                diff_origin[1] = molecule->coord_y[iatom];
                diff_origin[2] = molecule->coord_z[iatom];

                double **grad = alloc_gradients(bra, ket);

                /*
                 * numerical gradients
                 */
                double t0 = abs_time();
                libgrpp_overlap_integrals_num_grad(bra, ket, diff_origin, grad);
                time_num += abs_time() - t0;

                printf("numerical gradients:\n");
                print_gradients(bra, ket, grad);

                /*
                 * analytic gradients
                 */
                t0 = abs_time();
                libgrpp_overlap_integrals_analytic_grad(bra, ket, diff_origin, grad);
                time_ana += abs_time() - t0;

                printf("analytic gradients:\n");
                print_gradients(bra, ket, grad);

                dealloc_gradients(grad);
            }

            joffset += ket_dim;
        }

        ioffset += bra_dim;
    }

    printf("\n");
    printf("Time for numerical gradients: %.3f sec\n", time_num);
    printf("Time for analytic gradients: %.3f sec\n", time_ana);
    printf("\n");
}


/**
 * Gradients of the overlap integrals for a given shell pair
 * with respect to the point 'nuc_coord'
 * by the numerical differentiation.
 */
void libgrpp_overlap_integrals_num_grad(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        double *nuc_coord,
        double **grad
)
{
    const double h = 1e-6;

    int cart_size_A = shell_A->cart_size;
    int cart_size_B = shell_B->cart_size;
    int buf_size = cart_size_A * cart_size_B;

    /*
     * initializations: set gradients to zero
     */
    for (int icoord = 0; icoord < 3; icoord++) {
        memset(grad[icoord], 0, sizeof(double) * buf_size);
    }

    /*
     * d<A|A>/d... = 0
     */
    if (points_are_equal(shell_A->origin, shell_B->origin)) {
        return;
    }

    /*
     * d<A|B>/dC = 0
     */
    if (!points_are_equal(shell_A->origin, nuc_coord) &&
        !points_are_equal(shell_B->origin, nuc_coord)) {
        return;
    }

    double *buf_plus = (double *) calloc(buf_size, sizeof(double));
    double *buf_minus = (double *) calloc(buf_size, sizeof(double));

    // loop over X, Y, Z
    for (int icoord = 0; icoord < 3; icoord++) {

        // < df/dA | g >
        if (points_are_equal(shell_A->origin, nuc_coord)) {
            // A + h
            shell_A->origin[icoord] += h;
            libgrpp_overlap_integrals(shell_A, shell_B, buf_plus);
            shell_A->origin[icoord] -= h;
            // A - h
            shell_A->origin[icoord] -= h;
            libgrpp_overlap_integrals(shell_A, shell_B, buf_minus);
            shell_A->origin[icoord] += h;
            // finite difference formula
            for (int i = 0; i < buf_size; i++) {
                grad[icoord][i] += (buf_plus[i] - buf_minus[i]) / (2 * h);
            }
        }

        // < f | dg/dA >
        if (points_are_equal(shell_B->origin, nuc_coord)) {
            // A + h
            shell_B->origin[icoord] += h;
            libgrpp_overlap_integrals(shell_A, shell_B, buf_plus);
            shell_B->origin[icoord] -= h;
            // A - h
            shell_B->origin[icoord] -= h;
            libgrpp_overlap_integrals(shell_A, shell_B, buf_minus);
            shell_B->origin[icoord] += h;
            // finite difference formula
            for (int i = 0; i < buf_size; i++) {
                grad[icoord][i] += (buf_plus[i] - buf_minus[i]) / (2 * h);
            }
        }
    }

    free(buf_plus);
    free(buf_minus);
}


/**
 * Analytic calculation of gradients of the overlap integrals for a given shell pair
 * with respect to the point 'nuc_coord'.
 */
void libgrpp_overlap_integrals_analytic_grad(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        double *nuc_coord,
        double **grad
)
{
    int cart_size_A = shell_A->cart_size;
    int cart_size_B = shell_B->cart_size;
    int buf_size = cart_size_A * cart_size_B;

    /*
     * initializations: set gradients to zero
     */
    for (int icoord = 0; icoord < 3; icoord++) {
        memset(grad[icoord], 0, sizeof(double) * buf_size);
    }

    /*
     * integrals are zero:
     * (1) for 1-center integrals <A|A> (due to the translational invariance)
     * (2) d<A|B> / dC = 0 (integral is constant for the given 'nuc_coord')
     */
    if (points_are_equal(shell_A->origin, shell_B->origin)) {
        return;
    }

    /*
     * construct gradients:
     * d<A|B>/dA = + < df/dA | B >
     * d<A|B>/dB = - < df/dA | B >
     *
     * note that due to the property of translational invariance,
     * d<A|B>/dB = - d<A|B>/dA
     */
    if (points_are_equal(shell_A->origin, nuc_coord)) {
        overlap_gradient_diff_bra_contribution(shell_A, shell_B, grad, +1.0);
    }
    if (points_are_equal(shell_B->origin, nuc_coord)) {
        overlap_gradient_diff_bra_contribution(shell_A, shell_B, grad, -1.0);
    }
}


/**
 * Calculates contribution to gradients arising from the < df/dA | g > term:
 *
 * grad += factor * < df/dA | g >
 *
 * (bra basis function is differentiated).
 */
void overlap_gradient_diff_bra_contribution(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        double **grad,
        double factor
)
{
    /*
     * calculate overlap integrals < df/dA | B >
     */
    double *overlap_down = NULL;
    double *overlap_up = NULL;
    int cart_size_down = 0;
    int cart_size_up = 0;

    overlap_gradient_diff_bra_overlap_integrals(shell_A, shell_B, &overlap_down, &overlap_up, &cart_size_down,
                                                &cart_size_up);

    /*
     * construct contributions to gradients:
     * d<A|B>/dA += < df/dA | B >
     */
    for (int icoord = 0; icoord < 3; icoord++) {
        for (int i = 0; i < shell_A->cart_size; i++) {
            for (int j = 0; j < shell_B->cart_size; j++) {

                int *bra_nlm = shell_A->cart_list + 3 * i;
                int *ket_nlm = shell_B->cart_list + 3 * j;
                int index = i * shell_B->cart_size + j;

                /*
                 * contribution from the L-1 gaussian
                 */
                if (shell_A->L > 0) {
                    bra_nlm[icoord] -= 1;
                    int bra_index = nlm_to_linear(bra_nlm);
                    int ket_index = nlm_to_linear(ket_nlm);
                    bra_nlm[icoord] += 1;

                    grad[icoord][index] -=
                            factor * bra_nlm[icoord] * overlap_down[shell_B->cart_size * bra_index + ket_index];
                }

                /*
                 * contribution from the L+1 gaussian
                 */
                bra_nlm[icoord] += 1;
                int bra_index = nlm_to_linear(bra_nlm);
                int ket_index = nlm_to_linear(ket_nlm);
                bra_nlm[icoord] -= 1;

                grad[icoord][index] += factor * overlap_up[shell_B->cart_size * bra_index + ket_index];
            }
        }
    }

    if (overlap_down) {
        free(overlap_down);
    }
    free(overlap_up);
}


/**
 * To assemble the contribution < df/dA | g > to gradients, one have to differentiate
 * Gaussian function. Such a differentiation yields two Gaussians with angular momenta
 * L-1 ("down") and L+1 ("up"):
 * dG/dA -> G(L-1) and G(L+1)
 *
 * This function constructs overlap matrices with these "downgraded" and "upgraded"
 * Gaussian functions:
 * < G(L-1) | G' > and < G(L+1) | G' >
 *
 */
void overlap_gradient_diff_bra_overlap_integrals(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        double **overlap_down,
        double **overlap_up,
        int *cart_size_down,
        int *cart_size_up
)
{
    /*
     * differentiation of contracted Gaussian functions
     */
    libgrpp_shell_t *shell_A_down = NULL;
    libgrpp_shell_t *shell_A_up = NULL;
    differentiate_shell(shell_A, &shell_A_down, &shell_A_up);

    *cart_size_down = 0;
    if (shell_A_down != NULL) {
        *cart_size_down = shell_A_down->cart_size;
    }
    *cart_size_up = shell_A_up->cart_size;

    /*
     * overlap matrix:
     * < L-1 | L>
     */
    if (shell_A_down != NULL) {
        *overlap_down = (double *) calloc(shell_A_down->cart_size * shell_B->cart_size, sizeof(double));
        libgrpp_overlap_integrals(shell_A_down, shell_B, *overlap_down);
    }
    else {
        *overlap_down = NULL;
    }

    /*
     * overlap matrix:
     * < L+1 | L>
     */
    *overlap_up = (double *) calloc(shell_A_up->cart_size * shell_B->cart_size, sizeof(double));
    libgrpp_overlap_integrals(shell_A_up, shell_B, *overlap_up);

    /*
     * clean up
     */
    if (shell_A_down) {
        libgrpp_delete_shell(shell_A_down);
    }
    libgrpp_delete_shell(shell_A_up);
}


/**
 * Performs differentiation of a contracted Gaussian.
 *
 * Note that the "2 alpha" factors are absorbed into coefficients, while the 'n' factor is not.
 * The latter must be accounted for explicitly at the stage of gradient construction.
 * For more details, see:
 * T. Helgaker, P. Jorgensen, J. Olsen, Molecular Electronic-Structure Theory,
 * John Wiley & Sons Ltd, 2000.
 * Chapter 9.2.2, "Recurrence relations for Cartesian Gaussians"
 *
 */
void differentiate_shell(libgrpp_shell_t *shell, libgrpp_shell_t **shell_minus, libgrpp_shell_t **shell_plus)
{
    // downwards
    if (shell->L > 0) {
        *shell_minus = libgrpp_new_shell(shell->origin, shell->L - 1, shell->num_primitives, shell->coeffs,
                                         shell->alpha);

        for (int i = 0; i < shell->num_primitives; i++) {

            double alpha = shell->alpha[i];
            double L = shell->L;

            (*shell_minus)->coeffs[i] *= norm_factor(alpha, L) / norm_factor(alpha, L - 1);
        }
    }
    else {
        *shell_minus = NULL;
    }

    // upwards
    *shell_plus = libgrpp_new_shell(shell->origin, shell->L + 1, shell->num_primitives, shell->coeffs, shell->alpha);
    for (int i = 0; i < shell->num_primitives; i++) {

        double alpha = shell->alpha[i];
        double L = shell->L;

        (*shell_plus)->coeffs[i] *= 2.0 * alpha * norm_factor(alpha, L) / norm_factor(alpha, L + 1);
    }
}


/**
 * Checks if two 3d points coincide with each other.
 */
int points_are_equal(double *a, double *b)
{
    double const thresh = 1e-12;

    if (fabs(a[0] - b[0]) < thresh &&
        fabs(a[1] - b[1]) < thresh &&
        fabs(a[2] - b[2]) < thresh) {
        return 1;
    }

    return 0;
}


/**
 * calculates sequential ("linear") index of the (n,l,m) primitive in the cartesian shell
 */
int nlm_to_linear(int *nlm)
{
    int n = nlm[0];
    int l = nlm[1];
    int m = nlm[2];

    int L = n + l + m;
    int cart_size = (L + 1) * (L + 2) / 2;
    int *cart_list = libgrpp_generate_shell_cartesians(L);

    int index = 0;
    for (index = 0; index < cart_size; index++) {
        if (cart_list[3 * index + 0] == n &&
            cart_list[3 * index + 1] == l &&
            cart_list[3 * index + 2] == m) {
            break;
        }
    }

    free(cart_list);

    return index;
}


/**
 * Calculates normalization factor for the primitive Gaussian
 * with the exponential parameter 'alpha' and angular momentum L.
 */
double norm_factor(double alpha, int L)
{
    return pow(2 * alpha / M_PI, 0.75) * pow(4 * alpha, 0.5 * L);
}


/**
 * Prints beautiful table of gradients for a given shell pair.
 */
void print_gradients(libgrpp_shell_t *bra, libgrpp_shell_t *ket, double **grad)
{
    for (int i = 0; i < bra->cart_size; i++) {
        for (int j = 0; j < ket->cart_size; j++) {
            int nA = bra->cart_list[3 * i + 0];
            int lA = bra->cart_list[3 * i + 1];
            int mA = bra->cart_list[3 * i + 2];
            int nB = ket->cart_list[3 * j + 0];
            int lB = ket->cart_list[3 * j + 1];
            int mB = ket->cart_list[3 * j + 2];

            int index = i * ket->cart_size + j;

            printf("  %1d%1d%1d - %1d%1d%1d  %24.16f%24.16f%24.16f\n",
                   nA, lA, mA, nB, lB, mB,
                   grad[0][index], grad[1][index], grad[2][index]);
        }
    }
}


/**
 * Allocates memory for gradients for a given shell pair.
 */
double **alloc_gradients(libgrpp_shell_t *bra, libgrpp_shell_t *ket)
{
    size_t size = bra->cart_size * ket->cart_size;

    double **grad = (double **) calloc(3, sizeof(double *));
    for (int icoord = 0; icoord < 3; icoord++) {
        grad[icoord] = (double *) calloc(size, sizeof(double));
    }

    return grad;
}


/**
 * Deallocates arrays containing gradients of AO integrals.
 */
void dealloc_gradients(double **grad)
{
    free(grad[0]);
    free(grad[1]);
    free(grad[2]);
    free(grad);
}
