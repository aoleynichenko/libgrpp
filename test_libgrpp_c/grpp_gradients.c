/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2022 Alexander Oleynichenko
 */

/*
 * This example shown how to calculate numerical and analytic gradients
 * of GRPP integrals using the LIBGRPP library.
 *
 * The property of translational invariance is used to derive analytic expressions.
 * No new types of integrals arise.
 */

#include "grpp_gradients.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "abs_time.h"
#include "rpp.h"
#include "shell_list.h"


void libgrpp_grpp_integrals_num_grad(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        grpp_t *grpp_operator,
        double *grpp_origin,
        double *nuc_coord,
        double **grad_arep,
        double **grad_so_x,
        double **grad_so_y,
        double **grad_so_z
);


void libgrpp_grpp_integrals_analytic_grad(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        grpp_t *grpp_operator,
        double *grpp_origin,
        double *nuc_coord,
        double **grad_arep,
        double **grad_so_x,
        double **grad_so_y,
        double **grad_so_z
);


void evaluate_grpp_integrals_shell_pair(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        grpp_t *grpp_operator,
        double *grpp_origin,
        double *arep_matrix,
        double *so_x_matrix,
        double *so_y_matrix,
        double *so_z_matrix
);


void grpp_gradient_diff_bra_contribution(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        grpp_t *grpp_operator,
        double *grpp_origin,
        double **grad_arep,
        double **grad_so_x,
        double **grad_so_y,
        double **grad_so_z,
        double factor
);


void grpp_gradient_diff_bra_grpp_integrals(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        grpp_t *grpp_operator,
        double *grpp_origin,
        double **arep_matrix_down,
        double **so_x_matrix_down,
        double **so_y_matrix_down,
        double **so_z_matrix_down,
        double **arep_matrix_up,
        double **so_x_matrix_up,
        double **so_y_matrix_up,
        double **so_z_matrix_up,
        int *cart_size_down,
        int *cart_size_up
);


void grpp_gradient_diff_ket_contribution(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        grpp_t *grpp_operator,
        double *grpp_origin,
        double **grad_arep,
        double **grad_so_x,
        double **grad_so_y,
        double **grad_so_z,
        double factor
);


void grpp_gradient_diff_ket_grpp_integrals(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        grpp_t *grpp_operator,
        double *grpp_origin,
        double **arep_matrix_down,
        double **so_x_matrix_down,
        double **so_y_matrix_down,
        double **so_z_matrix_down,
        double **arep_matrix_up,
        double **so_x_matrix_up,
        double **so_y_matrix_up,
        double **so_z_matrix_up,
        int *cart_size_down,
        int *cart_size_up
);

void grpp_gradient_contribution(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        grpp_t *grpp_operator,
        double *grpp_origin,
        double **grad_arep,
        double **grad_so_x,
        double **grad_so_y,
        double **grad_so_z,
        int diff_bra,
        double factor
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
void evaluate_grpp_integrals_gradients(
        int num_shells,
        libgrpp_shell_t **shell_list,
        molecule_t *molecule,
        grpp_t **grpp_list
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

                for (int irpp = 0; irpp < molecule->n_atoms; irpp++) {
                    int z = molecule->charges[irpp];
                    grpp_t *grpp_operator = grpp_list[z];
                    if (grpp_operator == NULL) {
                        continue;
                    }

                    double grpp_origin[3];
                    grpp_origin[0] = molecule->coord_x[irpp];
                    grpp_origin[1] = molecule->coord_y[irpp];
                    grpp_origin[2] = molecule->coord_z[irpp];

                    printf("grpp origin %6.3f%6.3f%6.3f:\n",
                           molecule->coord_x[irpp], molecule->coord_y[irpp], molecule->coord_z[irpp]);

                    double **grad_arep = alloc_gradients(bra, ket);
                    double **grad_so_x = alloc_gradients(bra, ket);
                    double **grad_so_y = alloc_gradients(bra, ket);
                    double **grad_so_z = alloc_gradients(bra, ket);

                    /*
                     * numerical gradients
                     */
                    double t0 = abs_time();
                    libgrpp_grpp_integrals_num_grad(bra, ket, grpp_operator, grpp_origin, diff_origin, grad_arep,
                                                    grad_so_x, grad_so_y, grad_so_z);
                    time_num += abs_time() - t0;

                    printf("numerical gradients:\n");
                    print_gradients(bra, ket, grad_so_y);

                    /*
                     * analytic gradients
                     */
                    t0 = abs_time();
                    libgrpp_grpp_integrals_analytic_grad(bra, ket, grpp_operator, grpp_origin, diff_origin, grad_arep,
                                                         grad_so_x, grad_so_y, grad_so_z);
                    time_ana += abs_time() - t0;

                    printf("analytic gradients:\n");
                    print_gradients(bra, ket, grad_so_y);

                    dealloc_gradients(grad_arep);
                    dealloc_gradients(grad_so_x);
                    dealloc_gradients(grad_so_y);
                    dealloc_gradients(grad_so_z);

                }
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


void libgrpp_grpp_integrals_num_grad(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        grpp_t *grpp_operator,
        double *grpp_origin,
        double *nuc_coord,
        double **grad_arep,
        double **grad_so_x,
        double **grad_so_y,
        double **grad_so_z
)
{
    const double h = 1e-2;

    int cart_size_A = shell_A->cart_size;
    int cart_size_B = shell_B->cart_size;
    int buf_size = cart_size_A * cart_size_B;

    /*
     * initialization: set all gradients to zero
     */
    for (int icoord = 0; icoord < 3; icoord++) {
        memset(grad_arep[icoord], 0, sizeof(double) * buf_size);
        memset(grad_so_x[icoord], 0, sizeof(double) * buf_size);
        memset(grad_so_y[icoord], 0, sizeof(double) * buf_size);
        memset(grad_so_z[icoord], 0, sizeof(double) * buf_size);
    }

    /*
     * d<AAA>/d... = 0
     */
    if (points_are_equal(shell_A->origin, grpp_origin) &&
        points_are_equal(shell_B->origin, grpp_origin)) {
        return;
    }

    /*
     * d<ACB>/dD = 0
     */
    if (!points_are_equal(shell_A->origin, nuc_coord) &&
        !points_are_equal(shell_B->origin, nuc_coord) &&
        !points_are_equal(grpp_origin, nuc_coord)) {
        return;
    }

    double *buf_arep_plus = (double *) calloc(buf_size, sizeof(double));
    double *buf_so_x_plus = (double *) calloc(buf_size, sizeof(double));
    double *buf_so_y_plus = (double *) calloc(buf_size, sizeof(double));
    double *buf_so_z_plus = (double *) calloc(buf_size, sizeof(double));
    double *buf_arep_minus = (double *) calloc(buf_size, sizeof(double));
    double *buf_so_x_minus = (double *) calloc(buf_size, sizeof(double));
    double *buf_so_y_minus = (double *) calloc(buf_size, sizeof(double));
    double *buf_so_z_minus = (double *) calloc(buf_size, sizeof(double));

    int A_eq_D = points_are_equal(shell_A->origin, nuc_coord);
    int B_eq_D = points_are_equal(shell_B->origin, nuc_coord);
    int C_eq_D = points_are_equal(grpp_origin, nuc_coord);

    // loop over X, Y, Z
    for (int icoord = 0; icoord < 3; icoord++) {


        if (A_eq_D) {
            printf("A += h\n");
            shell_A->origin[icoord] += h;
        }
        if (B_eq_D && shell_A != shell_B) {
            printf("B += h\n");
            shell_B->origin[icoord] += h;
        }
        if (C_eq_D) {
            printf("C += h\n");
            grpp_origin[icoord] += h;
        }

        printf("A = %20.12f%20.12f%20.12f\n", shell_A->origin[0], shell_A->origin[1], shell_A->origin[2]);
        printf("B = %20.12f%20.12f%20.12f\n", shell_B->origin[0], shell_B->origin[1], shell_B->origin[2]);
        printf("C = %20.12f%20.12f%20.12f\n", grpp_origin[0], grpp_origin[1], grpp_origin[2]);

        evaluate_grpp_integrals_shell_pair(
                shell_A, shell_B, grpp_operator, grpp_origin,
                buf_arep_plus, buf_so_x_plus, buf_so_y_plus, buf_so_z_plus
        );
        printf("+h: %20.12f\n", buf_so_y_plus[0]);

        if (A_eq_D) {
            shell_A->origin[icoord] -= 2 * h;
        }
        if (B_eq_D && shell_A != shell_B) {
            shell_B->origin[icoord] -= 2 * h;
        }
        if (C_eq_D) {
            grpp_origin[icoord] -= 2 * h;
        }

        printf("A = %20.12f%20.12f%20.12f\n", shell_A->origin[0], shell_A->origin[1], shell_A->origin[2]);
        printf("B = %20.12f%20.12f%20.12f\n", shell_B->origin[0], shell_B->origin[1], shell_B->origin[2]);
        printf("C = %20.12f%20.12f%20.12f\n", grpp_origin[0], grpp_origin[1], grpp_origin[2]);

        evaluate_grpp_integrals_shell_pair(
                shell_A, shell_B, grpp_operator, grpp_origin,
                buf_arep_minus, buf_so_x_minus, buf_so_y_minus, buf_so_z_minus
        );
        printf("-h: %20.12f\n", buf_so_y_minus[0]);

        if (A_eq_D) {
            shell_A->origin[icoord] += h;
        }
        if (B_eq_D && shell_A != shell_B) {
            shell_B->origin[icoord] += h;
        }
        if (C_eq_D) {
            grpp_origin[icoord] += h;
        }

        // finite difference formula
        for (int i = 0; i < buf_size; i++) {
            grad_arep[icoord][i] += (buf_arep_plus[i] - buf_arep_minus[i]) / (2 * h);
            grad_so_x[icoord][i] += (buf_so_x_plus[i] - buf_so_x_minus[i]) / (2 * h);
            grad_so_y[icoord][i] += (buf_so_y_plus[i] - buf_so_y_minus[i]) / (2 * h);
            grad_so_z[icoord][i] += (buf_so_z_plus[i] - buf_so_z_minus[i]) / (2 * h);
        }
    }

    free(buf_arep_plus);
    free(buf_so_x_plus);
    free(buf_so_y_plus);
    free(buf_so_z_plus);
    free(buf_arep_minus);
    free(buf_so_x_minus);
    free(buf_so_y_minus);
    free(buf_so_z_minus);
}


void libgrpp_grpp_integrals_analytic_grad(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        grpp_t *grpp_operator,
        double *grpp_origin,
        double *nuc_coord,
        double **grad_arep,
        double **grad_so_x,
        double **grad_so_y,
        double **grad_so_z
)
{
    int cart_size_A = shell_A->cart_size;
    int cart_size_B = shell_B->cart_size;
    int buf_size = cart_size_A * cart_size_B;

    /*
     * initialization: set all gradients to zero
     */
    for (int icoord = 0; icoord < 3; icoord++) {
        memset(grad_arep[icoord], 0, sizeof(double) * buf_size);
        memset(grad_so_x[icoord], 0, sizeof(double) * buf_size);
        memset(grad_so_y[icoord], 0, sizeof(double) * buf_size);
        memset(grad_so_z[icoord], 0, sizeof(double) * buf_size);
    }

    /*
     * d<AAA>/d... = 0
     */
    if (points_are_equal(shell_A->origin, grpp_origin) &&
        points_are_equal(shell_B->origin, grpp_origin)) {
        return;
    }

    /*
     * d<ACB>/dD = 0
     */
    if (!points_are_equal(shell_A->origin, nuc_coord) &&
        !points_are_equal(shell_B->origin, nuc_coord) &&
        !points_are_equal(grpp_origin, nuc_coord)) {
        return;
    }


    double *A = shell_A->origin;
    double *B = shell_B->origin;
    double *C = grpp_origin;
    double *D = nuc_coord;

    const int diff_bra = 1;
    const int diff_ket = 0;

    /*
     * Type ACB
     */
    if (!points_are_equal(A, C) && !points_are_equal(C, B) && !points_are_equal(A, B)) {
        printf("ACB\n");
        if (points_are_equal(A, D)) {
            grpp_gradient_contribution(shell_A, shell_B, grpp_operator, grpp_origin,
                                       grad_arep, grad_so_x, grad_so_y, grad_so_z, diff_bra, +1.0);
        }
        if (points_are_equal(B, D)) {
            grpp_gradient_contribution(shell_A, shell_B, grpp_operator, grpp_origin,
                                       grad_arep, grad_so_x, grad_so_y, grad_so_z, diff_ket, +1.0);
        }
        if (points_are_equal(C, D)) {
            grpp_gradient_contribution(shell_A, shell_B, grpp_operator, grpp_origin,
                                       grad_arep, grad_so_x, grad_so_y, grad_so_z, diff_bra, -1.0);
            grpp_gradient_contribution(shell_A, shell_B, grpp_operator, grpp_origin,
                                       grad_arep, grad_so_x, grad_so_y, grad_so_z, diff_ket, -1.0);
        }
    }

    /*
     * Type ACA
     */
    if (points_are_equal(A, B) && !points_are_equal(A, C)) {
        printf("ACA\n");
        if (points_are_equal(A, D)) {
            grpp_gradient_diff_bra_contribution(shell_A, shell_B, grpp_operator, grpp_origin,
                                                grad_arep, grad_so_x, grad_so_y, grad_so_z, +1.0);
            printf("grad = %20.12f%20.12f%20.12f\n", grad_arep[0][0], grad_arep[1][0], grad_arep[2][0]);
            grpp_gradient_diff_ket_contribution(shell_A, shell_B, grpp_operator, grpp_origin,
                                                grad_arep, grad_so_x, grad_so_y, grad_so_z, +1.0);
            printf("grad = %20.12f%20.12f%20.12f\n", grad_arep[0][0], grad_arep[1][0], grad_arep[2][0]);
            /*grpp_gradient_contribution(shell_A, shell_B, grpp_operator, grpp_origin,
                                       grad_arep, grad_so_x, grad_so_y, grad_so_z, diff_bra, +1.0);
            grpp_gradient_contribution(shell_A, shell_B, grpp_operator, grpp_origin,
                                       grad_arep, grad_so_x, grad_so_y, grad_so_z, diff_ket, +1.0);*/
        }
        else {
            grpp_gradient_diff_bra_contribution(shell_A, shell_B, grpp_operator, grpp_origin,
                                       grad_arep, grad_so_x, grad_so_y, grad_so_z, -1.0);
            grpp_gradient_diff_ket_contribution(shell_A, shell_B, grpp_operator, grpp_origin,
                                       grad_arep, grad_so_x, grad_so_y, grad_so_z, -1.0);
            /*
            grpp_gradient_contribution(shell_A, shell_B, grpp_operator, grpp_origin,
                                       grad_arep, grad_so_x, grad_so_y, grad_so_z, diff_bra, -1.0);
            grpp_gradient_contribution(shell_A, shell_B, grpp_operator, grpp_origin,
                                       grad_arep, grad_so_x, grad_so_y, grad_so_z, diff_ket, -1.0);
             */
        }
    }

    /*
     * Type ACC
     */
    if (!points_are_equal(A, C) && points_are_equal(C, B)) {
        printf("ACC\n");
        if (points_are_equal(A, D)) {
            grpp_gradient_contribution(shell_A, shell_B, grpp_operator, grpp_origin,
                                       grad_arep, grad_so_x, grad_so_y, grad_so_z, diff_bra, +1.0);
        }
        else {
            grpp_gradient_contribution(shell_A, shell_B, grpp_operator, grpp_origin,
                                       grad_arep, grad_so_x, grad_so_y, grad_so_z, diff_bra, -1.0);
        }
    }

    /*
     * Type CCB
     */
    if (points_are_equal(A, C) && !points_are_equal(C, B)) {
        printf("CCB\n");
        if (points_are_equal(B, D)) {
            grpp_gradient_contribution(shell_A, shell_B, grpp_operator, grpp_origin,
                                       grad_arep, grad_so_x, grad_so_y, grad_so_z, diff_ket, +1.0);
        }
        else {
            grpp_gradient_contribution(shell_A, shell_B, grpp_operator, grpp_origin,
                                       grad_arep, grad_so_x, grad_so_y, grad_so_z, diff_ket, -1.0);
        }
    }
}


/**
 * Calculates contribution to gradients arising from the < df/dA | g > term:
 *
 * grad += factor * < df/dA | g >
 *
 * (bra basis function is differentiated).
 */
void grpp_gradient_diff_bra_contribution(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        grpp_t *grpp_operator,
        double *grpp_origin,
        double **grad_arep,
        double **grad_so_x,
        double **grad_so_y,
        double **grad_so_z,
        double factor
)
{
    /*
     * calculate overlap integrals < df/dA | B >
     */
    double *arep_matrix_down = NULL;
    double *so_x_matrix_down = NULL;
    double *so_y_matrix_down = NULL;
    double *so_z_matrix_down = NULL;
    double *arep_matrix_up = NULL;
    double *so_x_matrix_up = NULL;
    double *so_y_matrix_up = NULL;
    double *so_z_matrix_up = NULL;
    int cart_size_down = 0;
    int cart_size_up = 0;

    grpp_gradient_diff_bra_grpp_integrals(
            shell_A, shell_B, grpp_operator, grpp_origin,
            &arep_matrix_down, &so_x_matrix_down, &so_y_matrix_down, &so_z_matrix_down,
            &arep_matrix_up, &so_x_matrix_up, &so_y_matrix_up, &so_z_matrix_up,
            &cart_size_down, &cart_size_up
    );

    printf("bra contrib:\n");
    printf("%20.12f\n", so_y_matrix_up[0]);
    printf("%20.12f\n", so_y_matrix_up[1]);
    printf("%20.12f\n", so_y_matrix_up[2]);

    /*
     * construct contributions to gradients:
     * d<A|B>/dA += < df/dA | B >
     */
    for (int icoord = 0; icoord < 3; icoord++) {
        for (int i = 0; i < shell_A->cart_size; i++) {
            for (int j = 0; j < shell_B->cart_size; j++) {

                int bra_nlm[3];
                bra_nlm[0] = shell_A->cart_list[3 * i + 0];
                bra_nlm[1] = shell_A->cart_list[3 * i + 1];
                bra_nlm[2] = shell_A->cart_list[3 * i + 2];

                int ket_nlm[3];
                ket_nlm[0] = shell_B->cart_list[3 * j + 0];
                ket_nlm[1] = shell_B->cart_list[3 * j + 1];
                ket_nlm[2] = shell_B->cart_list[3 * j + 2];

                int index = i * shell_B->cart_size + j;

                /*
                 * contribution from the L-1 gaussian
                 */
                if (shell_A->L > 0) {
                    bra_nlm[icoord] -= 1;
                    int bra_index = nlm_to_linear(bra_nlm);
                    int ket_index = nlm_to_linear(ket_nlm);
                    bra_nlm[icoord] += 1;

                    grad_arep[icoord][index] -=
                            factor * bra_nlm[icoord] * arep_matrix_down[shell_B->cart_size * bra_index + ket_index];
                    grad_so_x[icoord][index] -=
                            factor * bra_nlm[icoord] * so_x_matrix_down[shell_B->cart_size * bra_index + ket_index];
                    grad_so_y[icoord][index] -=
                            factor * bra_nlm[icoord] * so_y_matrix_down[shell_B->cart_size * bra_index + ket_index];
                    grad_so_z[icoord][index] -=
                            factor * bra_nlm[icoord] * so_z_matrix_down[shell_B->cart_size * bra_index + ket_index];
                }

                /*
                 * contribution from the L+1 gaussian
                 */
                printf("bra_nlm = %d%d%d\n", bra_nlm[0], bra_nlm[1], bra_nlm[2]);
                bra_nlm[icoord] += 1;
                int bra_index = nlm_to_linear(bra_nlm);
                int ket_index = nlm_to_linear(ket_nlm);
                bra_nlm[icoord] -= 1;
                printf("bra_index = %d\n", bra_index);
                printf("ket_index = %d\n", ket_index);
                printf("index = %d\n", index);
                printf("index up = %d\n", shell_B->cart_size * bra_index + ket_index);
                printf("value = %20.12f\n", so_y_matrix_up[shell_B->cart_size * bra_index + ket_index]);

                grad_arep[icoord][index] += factor * arep_matrix_up[shell_B->cart_size * bra_index + ket_index];
                grad_so_x[icoord][index] += factor * so_x_matrix_up[shell_B->cart_size * bra_index + ket_index];
                grad_so_y[icoord][index] += factor * so_y_matrix_up[shell_B->cart_size * bra_index + ket_index];
                grad_so_z[icoord][index] += factor * so_z_matrix_up[shell_B->cart_size * bra_index + ket_index];
            }
        }
    }

    printf("grad value now:\n");
    printf("%20.12f\n", grad_so_y[0][0]);
    printf("%20.12f\n", grad_so_y[1][0]);
    printf("%20.12f\n", grad_so_y[2][0]);

    if (arep_matrix_down) {
        free(arep_matrix_down);
        free(so_x_matrix_down);
        free(so_y_matrix_down);
        free(so_z_matrix_down);
    }
    free(arep_matrix_up);
    free(so_x_matrix_up);
    free(so_y_matrix_up);
    free(so_z_matrix_up);
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
void grpp_gradient_diff_bra_grpp_integrals(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        grpp_t *grpp_operator,
        double *grpp_origin,
        double **arep_matrix_down,
        double **so_x_matrix_down,
        double **so_y_matrix_down,
        double **so_z_matrix_down,
        double **arep_matrix_up,
        double **so_x_matrix_up,
        double **so_y_matrix_up,
        double **so_z_matrix_up,
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
        size_t mat_size_down = shell_A_down->cart_size * shell_B->cart_size;
        *arep_matrix_down = (double *) calloc(mat_size_down, sizeof(double));
        *so_x_matrix_down = (double *) calloc(mat_size_down, sizeof(double));
        *so_y_matrix_down = (double *) calloc(mat_size_down, sizeof(double));
        *so_z_matrix_down = (double *) calloc(mat_size_down, sizeof(double));

        evaluate_grpp_integrals_shell_pair(
                shell_A_down, shell_B, grpp_operator, grpp_origin,
                *arep_matrix_down, *so_x_matrix_down, *so_y_matrix_down, *so_z_matrix_down
        );
    }
    else {
        *arep_matrix_down = NULL;
        *so_x_matrix_down = NULL;
        *so_y_matrix_down = NULL;
        *so_z_matrix_down = NULL;
    }

    /*
     * overlap matrix:
     * < L+1 | L>
     */
    size_t mat_size_up = shell_A_up->cart_size * shell_B->cart_size;
    *arep_matrix_up = (double *) calloc(mat_size_up, sizeof(double));
    *so_x_matrix_up = (double *) calloc(mat_size_up, sizeof(double));
    *so_y_matrix_up = (double *) calloc(mat_size_up, sizeof(double));
    *so_z_matrix_up = (double *) calloc(mat_size_up, sizeof(double));

    evaluate_grpp_integrals_shell_pair(
            shell_A_up, shell_B, grpp_operator, grpp_origin,
            *arep_matrix_up, *so_x_matrix_up, *so_y_matrix_up, *so_z_matrix_up
    );

    /*
     * clean up
     */
    if (shell_A_down) {
        libgrpp_delete_shell(shell_A_down);
    }
    libgrpp_delete_shell(shell_A_up);
}


/**
 * Calculates contribution to gradients arising from the < df/dA | g > term:
 *
 * grad += factor * < f | dg/dA >
 *
 * (bra basis function is differentiated).
 */
void grpp_gradient_diff_ket_contribution(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        grpp_t *grpp_operator,
        double *grpp_origin,
        double **grad_arep,
        double **grad_so_x,
        double **grad_so_y,
        double **grad_so_z,
        double factor
)
{
    /*
     * calculate overlap integrals < df/dA | B >
     */
    double *arep_matrix_down = NULL;
    double *so_x_matrix_down = NULL;
    double *so_y_matrix_down = NULL;
    double *so_z_matrix_down = NULL;
    double *arep_matrix_up = NULL;
    double *so_x_matrix_up = NULL;
    double *so_y_matrix_up = NULL;
    double *so_z_matrix_up = NULL;
    int cart_size_down = 0;
    int cart_size_up = 0;

    grpp_gradient_diff_ket_grpp_integrals(
            shell_A, shell_B, grpp_operator, grpp_origin,
            &arep_matrix_down, &so_x_matrix_down, &so_y_matrix_down, &so_z_matrix_down,
            &arep_matrix_up, &so_x_matrix_up, &so_y_matrix_up, &so_z_matrix_up,
            &cart_size_down, &cart_size_up
    );

    printf("ket contrib:\n");
    printf("%20.12f\n", so_y_matrix_up[0]);
    printf("%20.12f\n", so_y_matrix_up[1]);
    printf("%20.12f\n", so_y_matrix_up[2]);

    /*
     * construct contributions to gradients:
     * d<A|B>/dA += < df/dA | B >
     */
    for (int icoord = 0; icoord < 3; icoord++) {
        for (int i = 0; i < shell_A->cart_size; i++) {
            for (int j = 0; j < shell_B->cart_size; j++) {

                int bra_nlm[3];
                bra_nlm[0] = shell_A->cart_list[3 * i + 0];
                bra_nlm[1] = shell_A->cart_list[3 * i + 1];
                bra_nlm[2] = shell_A->cart_list[3 * i + 2];

                int ket_nlm[3];
                ket_nlm[0] = shell_B->cart_list[3 * j + 0];
                ket_nlm[1] = shell_B->cart_list[3 * j + 1];
                ket_nlm[2] = shell_B->cart_list[3 * j + 2];

                int index = i * shell_B->cart_size + j;

                /*
                 * contribution from the L-1 gaussian
                 */
                if (shell_B->L > 0) {
                    ket_nlm[icoord] -= 1;
                    int bra_index = nlm_to_linear(bra_nlm);
                    int ket_index = nlm_to_linear(ket_nlm);
                    ket_nlm[icoord] += 1;

                    grad_arep[icoord][index] -=
                            factor * ket_nlm[icoord] * arep_matrix_down[cart_size_down * bra_index + ket_index];
                    grad_so_x[icoord][index] -=
                            factor * ket_nlm[icoord] * so_x_matrix_down[cart_size_down * bra_index + ket_index];
                    grad_so_y[icoord][index] -=
                            factor * ket_nlm[icoord] * so_y_matrix_down[cart_size_down * bra_index + ket_index];
                    grad_so_z[icoord][index] -=
                            factor * ket_nlm[icoord] * so_z_matrix_down[cart_size_down * bra_index + ket_index];
                }

                /*
                 * contribution from the L+1 gaussian
                 */
                printf("ket_nlm = %d%d%d\n", ket_nlm[0], ket_nlm[1], ket_nlm[2]);

                ket_nlm[icoord] += 1;
                int bra_index = nlm_to_linear(bra_nlm);
                int ket_index = nlm_to_linear(ket_nlm);
                ket_nlm[icoord] -= 1;

                printf("bra_index = %d\n", bra_index);
                printf("ket_index = %d\n", ket_index);
                printf("index = %d\n", index);
                printf("index up = %d\n", cart_size_up * bra_index + ket_index);
                printf("value = %20.12f\n", so_y_matrix_up[cart_size_up * bra_index + ket_index]);

                grad_arep[icoord][index] += factor * arep_matrix_up[cart_size_up * bra_index + ket_index];
                grad_so_x[icoord][index] += factor * so_x_matrix_up[cart_size_up * bra_index + ket_index];
                grad_so_y[icoord][index] += factor * so_y_matrix_up[cart_size_up * bra_index + ket_index];
                grad_so_z[icoord][index] += factor * so_z_matrix_up[cart_size_up * bra_index + ket_index];
            }
        }
    }

    if (arep_matrix_down) {
        free(arep_matrix_down);
        free(so_x_matrix_down);
        free(so_y_matrix_down);
        free(so_z_matrix_down);
    }
    free(arep_matrix_up);
    free(so_x_matrix_up);
    free(so_y_matrix_up);
    free(so_z_matrix_up);
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
void grpp_gradient_diff_ket_grpp_integrals(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        grpp_t *grpp_operator,
        double *grpp_origin,
        double **arep_matrix_down,
        double **so_x_matrix_down,
        double **so_y_matrix_down,
        double **so_z_matrix_down,
        double **arep_matrix_up,
        double **so_x_matrix_up,
        double **so_y_matrix_up,
        double **so_z_matrix_up,
        int *cart_size_down,
        int *cart_size_up
)
{
    /*
     * differentiation of contracted Gaussian functions
     */
    libgrpp_shell_t *shell_B_down = NULL;
    libgrpp_shell_t *shell_B_up = NULL;
    differentiate_shell(shell_B, &shell_B_down, &shell_B_up);

    *cart_size_down = 0;
    if (shell_B_down != NULL) {
        *cart_size_down = shell_B_down->cart_size;
    }
    *cart_size_up = shell_B_up->cart_size;

    /*
     * overlap matrix:
     * < L-1 | L>
     */
    if (shell_B_down != NULL) {
        size_t mat_size_down = shell_A->cart_size * shell_B_down->cart_size;
        *arep_matrix_down = (double *) calloc(mat_size_down, sizeof(double));
        *so_x_matrix_down = (double *) calloc(mat_size_down, sizeof(double));
        *so_y_matrix_down = (double *) calloc(mat_size_down, sizeof(double));
        *so_z_matrix_down = (double *) calloc(mat_size_down, sizeof(double));

        evaluate_grpp_integrals_shell_pair(
                shell_A, shell_B_down, grpp_operator, grpp_origin,
                *arep_matrix_down, *so_x_matrix_down, *so_y_matrix_down, *so_z_matrix_down
        );
    }
    else {
        *arep_matrix_down = NULL;
        *so_x_matrix_down = NULL;
        *so_y_matrix_down = NULL;
        *so_z_matrix_down = NULL;
    }

    /*
     * overlap matrix:
     * < L+1 | L>
     */
    size_t mat_size_up = shell_A->cart_size * shell_B_up->cart_size;
    *arep_matrix_up = (double *) calloc(mat_size_up, sizeof(double));
    *so_x_matrix_up = (double *) calloc(mat_size_up, sizeof(double));
    *so_y_matrix_up = (double *) calloc(mat_size_up, sizeof(double));
    *so_z_matrix_up = (double *) calloc(mat_size_up, sizeof(double));

    evaluate_grpp_integrals_shell_pair(
            shell_A, shell_B_up, grpp_operator, grpp_origin,
            *arep_matrix_up, *so_x_matrix_up, *so_y_matrix_up, *so_z_matrix_up
    );

    /*
     * clean up
     */
    if (shell_B_down) {
        libgrpp_delete_shell(shell_B_down);
    }
    libgrpp_delete_shell(shell_B_up);
}


void grpp_gradient_diff_gaussian(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        grpp_t *grpp_operator,
        double *grpp_origin,
        double **arep_matrix_down,
        double **so_x_matrix_down,
        double **so_y_matrix_down,
        double **so_z_matrix_down,
        double **arep_matrix_up,
        double **so_x_matrix_up,
        double **so_y_matrix_up,
        double **so_z_matrix_up,
        int *cart_size_down,
        int *cart_size_up,
        int diff_bra
);


void grpp_gradient_contribution(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        grpp_t *grpp_operator,
        double *grpp_origin,
        double **grad_arep,
        double **grad_so_x,
        double **grad_so_y,
        double **grad_so_z,
        int diff_bra,
        double factor
)
{
    int diff_ket = 0;

    if (diff_bra == 0) {
        diff_bra = 0;
        diff_ket = 1;
    }
    else {
        diff_bra = 1;
        diff_ket = 0;
    }

    /*
     * calculate overlap integrals < df/dA | B >
     */
    double *arep_matrix_down = NULL;
    double *so_x_matrix_down = NULL;
    double *so_y_matrix_down = NULL;
    double *so_z_matrix_down = NULL;
    double *arep_matrix_up = NULL;
    double *so_x_matrix_up = NULL;
    double *so_y_matrix_up = NULL;
    double *so_z_matrix_up = NULL;
    int cart_size_down = 0;
    int cart_size_up = 0;

    grpp_gradient_diff_gaussian(
            shell_A, shell_B, grpp_operator, grpp_origin,
            &arep_matrix_down, &so_x_matrix_down, &so_y_matrix_down, &so_z_matrix_down,
            &arep_matrix_up, &so_x_matrix_up, &so_y_matrix_up, &so_z_matrix_up,
            &cart_size_down, &cart_size_up, diff_bra
    );

    /*
     * construct contributions to gradients:
     * d<A|U|B>/dA += < df/dA | U | B >
     */
    for (int icoord = 0; icoord < 3; icoord++) {
        for (int i = 0; i < shell_A->cart_size; i++) {
            for (int j = 0; j < shell_B->cart_size; j++) {

                int bra_nlm[3];
                bra_nlm[0] = shell_A->cart_list[3 * i + 0];
                bra_nlm[1] = shell_A->cart_list[3 * i + 1];
                bra_nlm[2] = shell_A->cart_list[3 * i + 2];

                int ket_nlm[3];
                ket_nlm[0] = shell_B->cart_list[3 * j + 0];
                ket_nlm[1] = shell_B->cart_list[3 * j + 1];
                ket_nlm[2] = shell_B->cart_list[3 * j + 2];

                int index = i * shell_B->cart_size + j;

                int *diff_nlm = diff_bra ? bra_nlm : ket_nlm;

                /*
                 * contribution from the L-1 gaussian
                 */
                if (cart_size_down > 0) {
                    diff_nlm[icoord] -= 1;
                    int bra_index = nlm_to_linear(bra_nlm);
                    int ket_index = nlm_to_linear(ket_nlm);
                    diff_nlm[icoord] += 1;

                    int n = diff_nlm[icoord];
                    int row_len = diff_bra ? shell_B->cart_size : cart_size_down;
                    int index_down = row_len * bra_index + ket_index;

                    grad_arep[icoord][index] -= factor * n * arep_matrix_down[index_down];
                    grad_so_x[icoord][index] -= factor * n * so_x_matrix_down[index_down];
                    grad_so_y[icoord][index] -= factor * n * so_y_matrix_down[index_down];
                    grad_so_z[icoord][index] -= factor * n * so_z_matrix_down[index_down];
                }

                /*
                 * contribution from the L+1 gaussian
                 */
                diff_nlm[icoord] += 1;
                int bra_index = nlm_to_linear(bra_nlm);
                int ket_index = nlm_to_linear(ket_nlm);
                diff_nlm[icoord] -= 1;

                int row_len = diff_bra ? shell_B->cart_size : cart_size_up;
                int index_up = row_len * bra_index + ket_index;

                grad_arep[icoord][index] += factor * arep_matrix_up[index_up];
                grad_so_x[icoord][index] += factor * so_x_matrix_up[index_up];
                grad_so_y[icoord][index] += factor * so_y_matrix_up[index_up];
                grad_so_z[icoord][index] += factor * so_z_matrix_up[index_up];
            }
        }
    }

    if (arep_matrix_down) {
        free(arep_matrix_down);
        free(so_x_matrix_down);
        free(so_y_matrix_down);
        free(so_z_matrix_down);
    }
    free(arep_matrix_up);
    free(so_x_matrix_up);
    free(so_y_matrix_up);
    free(so_z_matrix_up);
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
void grpp_gradient_diff_gaussian(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        grpp_t *grpp_operator,
        double *grpp_origin,
        double **arep_matrix_down,
        double **so_x_matrix_down,
        double **so_y_matrix_down,
        double **so_z_matrix_down,
        double **arep_matrix_up,
        double **so_x_matrix_up,
        double **so_y_matrix_up,
        double **so_z_matrix_up,
        int *cart_size_down,
        int *cart_size_up,
        int diff_bra
)
{
    int diff_ket = 0;

    if (diff_bra == 0) {
        diff_bra = 0;
        diff_ket = 1;
    }
    else {
        diff_bra = 1;
        diff_ket = 0;
    }

    /*
     * which shell should be differentiated, bra or ket
     */
    libgrpp_shell_t *const_shell = NULL;
    libgrpp_shell_t *diff_shell = NULL;
    if (diff_bra) {
        diff_shell = shell_A;
        const_shell = shell_B;
    }
    else {
        diff_shell = shell_B;
        const_shell = shell_A;
    }

    /*
     * differentiation of contracted Gaussian functions
     */
    libgrpp_shell_t *shell_down = NULL;
    libgrpp_shell_t *shell_up = NULL;
    differentiate_shell(diff_shell, &shell_down, &shell_up);

    *cart_size_down = 0;
    if (shell_down != NULL) {
        *cart_size_down = shell_down->cart_size;
    }
    *cart_size_up = shell_up->cart_size;

    /*
     * GRPP matrix:
     * < L-1 | U | L > or < L | U | L-1 >
     */
    if (shell_down != NULL) {
        size_t mat_size_down = const_shell->cart_size * shell_down->cart_size;
        *arep_matrix_down = (double *) calloc(mat_size_down, sizeof(double));
        *so_x_matrix_down = (double *) calloc(mat_size_down, sizeof(double));
        *so_y_matrix_down = (double *) calloc(mat_size_down, sizeof(double));
        *so_z_matrix_down = (double *) calloc(mat_size_down, sizeof(double));

        evaluate_grpp_integrals_shell_pair(
                diff_bra ? shell_down : shell_A,
                diff_bra ? shell_B : shell_down,
                grpp_operator, grpp_origin,
                *arep_matrix_down, *so_x_matrix_down, *so_y_matrix_down, *so_z_matrix_down
        );
    }
    else {
        *arep_matrix_down = NULL;
        *so_x_matrix_down = NULL;
        *so_y_matrix_down = NULL;
        *so_z_matrix_down = NULL;
    }

    /*
     * GRPP matrix:
     * < L+1 | U | L > or < L | U | L+1 >
     */
    size_t mat_size_up = const_shell->cart_size * shell_up->cart_size;
    *arep_matrix_up = (double *) calloc(mat_size_up, sizeof(double));
    *so_x_matrix_up = (double *) calloc(mat_size_up, sizeof(double));
    *so_y_matrix_up = (double *) calloc(mat_size_up, sizeof(double));
    *so_z_matrix_up = (double *) calloc(mat_size_up, sizeof(double));

    evaluate_grpp_integrals_shell_pair(
            diff_bra ? shell_up : shell_A,
            diff_bra ? shell_B : shell_up,
            grpp_operator, grpp_origin,
            *arep_matrix_up, *so_x_matrix_up, *so_y_matrix_up, *so_z_matrix_up
    );

    /*
     * clean up
     */
    if (shell_down) {
        libgrpp_delete_shell(shell_down);
    }
    libgrpp_delete_shell(shell_up);
}

