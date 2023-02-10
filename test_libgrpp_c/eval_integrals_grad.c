/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2023 Alexander Oleynichenko
 */

#include "eval_integrals_grad.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "abs_time.h"
#include "../libgrpp/libgrpp.h"
#include "molecule.h"
#include "shell_list.h"


#define MAX_BUF 10000

void add_block_to_matrix(int dim_1, int dim_2, double *matrix,
                         int block_dim_1, int block_dim_2, double *block,
                         int col_offset, int row_offset, double factor);

void print_gradients(libgrpp_shell_t *bra, libgrpp_shell_t *ket, double **grad);

double **alloc_gradients(libgrpp_shell_t *bra, libgrpp_shell_t *ket);

void dealloc_gradients(double **grad);


/**
 * This function evaluates gradients of overlap integrals with respect to coordinates
 * of all atoms in a given molecule.
 */
void evaluate_overlap_integrals_gradient(
        int num_shells,
        libgrpp_shell_t **shell_list,
        molecule_t *molecule,
        double **gradient
)
{
    int dim = calculate_basis_dim(shell_list, num_shells);

    int ioffset = 0;
    for (int ishell = 0; ishell < num_shells; ishell++) {

        libgrpp_shell_t *bra = shell_list[ishell];
        int bra_dim = libgrpp_get_shell_size(bra);

        int joffset = 0;
        for (int jshell = 0; jshell < num_shells; jshell++) {

            libgrpp_shell_t *ket = shell_list[jshell];
            int ket_dim = libgrpp_get_shell_size(ket);

            double t1 = abs_time();
            printf("overlap grad: ishell=%3d (L=%d)\tjshell=%3d (L=%d)\t", ishell + 1, bra->L, jshell + 1, ket->L);

            double **buf_grad = alloc_gradients(bra, ket);

            for (int iatom = 0; iatom < molecule->n_atoms; iatom++) {

                double point_3d[3];
                point_3d[0] = molecule->coord_x[iatom];
                point_3d[1] = molecule->coord_y[iatom];
                point_3d[2] = molecule->coord_z[iatom];

                libgrpp_overlap_integrals_gradient(bra, ket, point_3d, buf_grad);

                double *grad_x = gradient[3 * iatom + 0];
                double *grad_y = gradient[3 * iatom + 1];
                double *grad_z = gradient[3 * iatom + 2];

                add_block_to_matrix(dim, dim, grad_x, bra_dim, ket_dim, buf_grad[0], ioffset, joffset, 1.0);
                add_block_to_matrix(dim, dim, grad_y, bra_dim, ket_dim, buf_grad[1], ioffset, joffset, 1.0);
                add_block_to_matrix(dim, dim, grad_z, bra_dim, ket_dim, buf_grad[2], ioffset, joffset, 1.0);
            }

            double t2 = abs_time();
            printf("%10.3f sec\n", t2 - t1);

            dealloc_gradients(buf_grad);

            joffset += ket_dim;
        }

        ioffset += bra_dim;
    }
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

