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
 * This function evaluates gradients of grpp integrals with respect to coordinates
 * of all atoms in a given molecule.
 */
void evaluate_grpp_integrals_gradient(
        int num_shells,
        libgrpp_shell_t **shell_list,
        molecule_t *molecule,
        libgrpp_grpp_t **grpp_list,
        double **gradient_arep,
        double **gradient_so_x,
        double **gradient_so_y,
        double **gradient_so_z
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
            printf("grpp grad: ishell=%3d (L=%d)\tjshell=%3d (L=%d)\t", ishell + 1, bra->L, jshell + 1, ket->L);

            for (int iatom = 0; iatom < molecule->n_atoms; iatom++) {

                double point_3d[3];
                point_3d[0] = molecule->coord_x[iatom];
                point_3d[1] = molecule->coord_y[iatom];
                point_3d[2] = molecule->coord_z[iatom];

                for (int irpp = 0; irpp < molecule->n_atoms; irpp++) {
                    int z = molecule->charges[irpp];
                    libgrpp_grpp_t *grpp_operator = grpp_list[z];
                    if (grpp_operator == NULL) {
                        continue;
                    }

                    double grpp_origin[3];
                    grpp_origin[0] = molecule->coord_x[irpp];
                    grpp_origin[1] = molecule->coord_y[irpp];
                    grpp_origin[2] = molecule->coord_z[irpp];

                    double **buf_grad_arep = alloc_gradients(bra, ket);
                    double **buf_grad_so_x = alloc_gradients(bra, ket);
                    double **buf_grad_so_y = alloc_gradients(bra, ket);
                    double **buf_grad_so_z = alloc_gradients(bra, ket);

                    libgrpp_full_grpp_integrals_gradient(bra, ket, grpp_operator, grpp_origin, point_3d, buf_grad_arep,
                                                         buf_grad_so_x, buf_grad_so_y, buf_grad_so_z);

                    /*
                     * AREP
                     */
                    double *grad_x_arep = gradient_arep[3 * iatom + 0];
                    double *grad_y_arep = gradient_arep[3 * iatom + 1];
                    double *grad_z_arep = gradient_arep[3 * iatom + 2];

                    add_block_to_matrix(dim, dim, grad_x_arep, bra_dim, ket_dim, buf_grad_arep[0], ioffset, joffset,
                                        1.0);
                    add_block_to_matrix(dim, dim, grad_y_arep, bra_dim, ket_dim, buf_grad_arep[1], ioffset, joffset,
                                        1.0);
                    add_block_to_matrix(dim, dim, grad_z_arep, bra_dim, ket_dim, buf_grad_arep[2], ioffset, joffset,
                                        1.0);

                    /*
                     * SO-X
                     */
                    double *grad_x_so_x = gradient_so_x[3 * iatom + 0];
                    double *grad_y_so_x = gradient_so_x[3 * iatom + 1];
                    double *grad_z_so_x = gradient_so_x[3 * iatom + 2];

                    add_block_to_matrix(dim, dim, grad_x_so_x, bra_dim, ket_dim, buf_grad_so_x[0], ioffset, joffset,
                                        1.0);
                    add_block_to_matrix(dim, dim, grad_y_so_x, bra_dim, ket_dim, buf_grad_so_x[1], ioffset, joffset,
                                        1.0);
                    add_block_to_matrix(dim, dim, grad_z_so_x, bra_dim, ket_dim, buf_grad_so_x[2], ioffset, joffset,
                                        1.0);

                    /*
                     * SO-Y
                     */
                    double *grad_x_so_y = gradient_so_y[3 * iatom + 0];
                    double *grad_y_so_y = gradient_so_y[3 * iatom + 1];
                    double *grad_z_so_y = gradient_so_y[3 * iatom + 2];

                    add_block_to_matrix(dim, dim, grad_x_so_y, bra_dim, ket_dim, buf_grad_so_y[0], ioffset, joffset,
                                        1.0);
                    add_block_to_matrix(dim, dim, grad_y_so_y, bra_dim, ket_dim, buf_grad_so_y[1], ioffset, joffset,
                                        1.0);
                    add_block_to_matrix(dim, dim, grad_z_so_y, bra_dim, ket_dim, buf_grad_so_y[2], ioffset, joffset,
                                        1.0);

                    /*
                     * SO-Z
                     */
                    double *grad_x_so_z = gradient_so_z[3 * iatom + 0];
                    double *grad_y_so_z = gradient_so_z[3 * iatom + 1];
                    double *grad_z_so_z = gradient_so_z[3 * iatom + 2];

                    add_block_to_matrix(dim, dim, grad_x_so_z, bra_dim, ket_dim, buf_grad_so_z[0], ioffset, joffset,
                                        1.0);
                    add_block_to_matrix(dim, dim, grad_y_so_z, bra_dim, ket_dim, buf_grad_so_z[1], ioffset, joffset,
                                        1.0);
                    add_block_to_matrix(dim, dim, grad_z_so_z, bra_dim, ket_dim, buf_grad_so_z[2], ioffset, joffset,
                                        1.0);

                    dealloc_gradients(buf_grad_arep);
                    dealloc_gradients(buf_grad_so_x);
                    dealloc_gradients(buf_grad_so_y);
                    dealloc_gradients(buf_grad_so_z);
                }
            }

            double t2 = abs_time();
            printf("%10.3f sec\n", t2 - t1);

            joffset += ket_dim;
        }

        ioffset += bra_dim;
    }
}


/**
 * This function evaluates gradients of grpp integrals with respect to coordinates
 * of all atoms in a given molecule.
 */
void evaluate_type1_integrals_gradient(
        int num_shells,
        libgrpp_shell_t **shell_list,
        molecule_t *molecule,
        libgrpp_grpp_t **grpp_list,
        double **gradient_arep
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
            printf("type1 grad: ishell=%3d (L=%d)\tjshell=%3d (L=%d)\t", ishell + 1, bra->L, jshell + 1, ket->L);

            for (int iatom = 0; iatom < molecule->n_atoms; iatom++) {

                double point_3d[3];
                point_3d[0] = molecule->coord_x[iatom];
                point_3d[1] = molecule->coord_y[iatom];
                point_3d[2] = molecule->coord_z[iatom];

                for (int irpp = 0; irpp < molecule->n_atoms; irpp++) {
                    int z = molecule->charges[irpp];
                    libgrpp_grpp_t *grpp_operator = grpp_list[z];
                    if (grpp_operator == NULL) {
                        continue;
                    }

                    double grpp_origin[3];
                    grpp_origin[0] = molecule->coord_x[irpp];
                    grpp_origin[1] = molecule->coord_y[irpp];
                    grpp_origin[2] = molecule->coord_z[irpp];

                    double **buf_grad_arep = alloc_gradients(bra, ket);
                    libgrpp_type1_integrals_gradient(
                            bra, ket, grpp_origin, grpp_operator->U_L, point_3d, buf_grad_arep
                    );

                    /*
                     * AREP
                     */
                    double *grad_x_arep = gradient_arep[3 * iatom + 0];
                    double *grad_y_arep = gradient_arep[3 * iatom + 1];
                    double *grad_z_arep = gradient_arep[3 * iatom + 2];

                    add_block_to_matrix(dim, dim, grad_x_arep, bra_dim, ket_dim, buf_grad_arep[0], ioffset, joffset,
                                        1.0);
                    add_block_to_matrix(dim, dim, grad_y_arep, bra_dim, ket_dim, buf_grad_arep[1], ioffset, joffset,
                                        1.0);
                    add_block_to_matrix(dim, dim, grad_z_arep, bra_dim, ket_dim, buf_grad_arep[2], ioffset, joffset,
                                        1.0);

                    dealloc_gradients(buf_grad_arep);
                }
            }

            double t2 = abs_time();
            printf("%10.3f sec\n", t2 - t1);

            joffset += ket_dim;
        }

        ioffset += bra_dim;
    }
}


/**
 * This function evaluates gradients of grpp integrals with respect to coordinates
 * of all atoms in a given molecule.
 */
void evaluate_type2_integrals_gradient(
        int num_shells,
        libgrpp_shell_t **shell_list,
        molecule_t *molecule,
        libgrpp_grpp_t **grpp_list,
        double **gradient_arep
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
            printf("type2 grad: ishell=%3d (L=%d)\tjshell=%3d (L=%d)\t", ishell + 1, bra->L, jshell + 1, ket->L);

            for (int iatom = 0; iatom < molecule->n_atoms; iatom++) {

                double point_3d[3];
                point_3d[0] = molecule->coord_x[iatom];
                point_3d[1] = molecule->coord_y[iatom];
                point_3d[2] = molecule->coord_z[iatom];

                for (int irpp = 0; irpp < molecule->n_atoms; irpp++) {
                    int z = molecule->charges[irpp];
                    libgrpp_grpp_t *grpp_operator = grpp_list[z];
                    if (grpp_operator == NULL) {
                        continue;
                    }

                    double grpp_origin[3];
                    grpp_origin[0] = molecule->coord_x[irpp];
                    grpp_origin[1] = molecule->coord_y[irpp];
                    grpp_origin[2] = molecule->coord_z[irpp];

                    for (int ipot = 0; ipot < grpp_operator->n_arep; ipot++) {

                        double **buf_grad_arep = alloc_gradients(bra, ket);
                        libgrpp_type2_integrals_gradient(
                                bra, ket, grpp_origin, grpp_operator->U_arep[ipot], point_3d, buf_grad_arep
                        );

                        /*
                         * AREP
                         */
                        double *grad_x_arep = gradient_arep[3 * iatom + 0];
                        double *grad_y_arep = gradient_arep[3 * iatom + 1];
                        double *grad_z_arep = gradient_arep[3 * iatom + 2];

                        add_block_to_matrix(dim, dim, grad_x_arep, bra_dim, ket_dim, buf_grad_arep[0], ioffset, joffset,
                                            1.0);
                        add_block_to_matrix(dim, dim, grad_y_arep, bra_dim, ket_dim, buf_grad_arep[1], ioffset, joffset,
                                            1.0);
                        add_block_to_matrix(dim, dim, grad_z_arep, bra_dim, ket_dim, buf_grad_arep[2], ioffset, joffset,
                                            1.0);

                        dealloc_gradients(buf_grad_arep);
                    }
                }
            }

            double t2 = abs_time();
            printf("%10.3f sec\n", t2 - t1);

            joffset += ket_dim;
        }

        ioffset += bra_dim;
    }
}


void evaluate_spin_orbit_integrals_gradient(
        int num_shells,
        libgrpp_shell_t **shell_list,
        molecule_t *molecule,
        libgrpp_grpp_t **grpp_list,
        double **gradient_so_x,
        double **gradient_so_y,
        double **gradient_so_z
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
            printf("spin-orbit grad: ishell=%3d (L=%d)\tjshell=%3d (L=%d)\t", ishell + 1, bra->L, jshell + 1, ket->L);

            for (int iatom = 0; iatom < molecule->n_atoms; iatom++) {

                double point_3d[3];
                point_3d[0] = molecule->coord_x[iatom];
                point_3d[1] = molecule->coord_y[iatom];
                point_3d[2] = molecule->coord_z[iatom];

                for (int irpp = 0; irpp < molecule->n_atoms; irpp++) {
                    int z = molecule->charges[irpp];
                    libgrpp_grpp_t *grpp_operator = grpp_list[z];
                    if (grpp_operator == NULL) {
                        continue;
                    }

                    double grpp_origin[3];
                    grpp_origin[0] = molecule->coord_x[irpp];
                    grpp_origin[1] = molecule->coord_y[irpp];
                    grpp_origin[2] = molecule->coord_z[irpp];

                    for (int ipot = 1; ipot < grpp_operator->n_esop; ipot++) {

                        double **buf_grad_so_x = alloc_gradients(bra, ket);
                        double **buf_grad_so_y = alloc_gradients(bra, ket);
                        double **buf_grad_so_z = alloc_gradients(bra, ket);

                        libgrpp_spin_orbit_integrals_gradient(
                                bra, ket, grpp_origin, grpp_operator->U_esop[ipot], point_3d,
                                buf_grad_so_x, buf_grad_so_y, buf_grad_so_z
                        );

                        /*
                         * SO-X
                         */
                        double *grad_x_so_x = gradient_so_x[3 * iatom + 0];
                        double *grad_y_so_x = gradient_so_x[3 * iatom + 1];
                        double *grad_z_so_x = gradient_so_x[3 * iatom + 2];

                        add_block_to_matrix(dim, dim, grad_x_so_x, bra_dim, ket_dim, buf_grad_so_x[0], ioffset, joffset,
                                            1.0);
                        add_block_to_matrix(dim, dim, grad_y_so_x, bra_dim, ket_dim, buf_grad_so_x[1], ioffset, joffset,
                                            1.0);
                        add_block_to_matrix(dim, dim, grad_z_so_x, bra_dim, ket_dim, buf_grad_so_x[2], ioffset, joffset,
                                            1.0);

                        /*
                         * SO-Y
                         */
                        double *grad_x_so_y = gradient_so_y[3 * iatom + 0];
                        double *grad_y_so_y = gradient_so_y[3 * iatom + 1];
                        double *grad_z_so_y = gradient_so_y[3 * iatom + 2];

                        add_block_to_matrix(dim, dim, grad_x_so_y, bra_dim, ket_dim, buf_grad_so_y[0], ioffset, joffset,
                                            1.0);
                        add_block_to_matrix(dim, dim, grad_y_so_y, bra_dim, ket_dim, buf_grad_so_y[1], ioffset, joffset,
                                            1.0);
                        add_block_to_matrix(dim, dim, grad_z_so_y, bra_dim, ket_dim, buf_grad_so_y[2], ioffset, joffset,
                                            1.0);

                        /*
                         * SO-Z
                         */
                        double *grad_x_so_z = gradient_so_z[3 * iatom + 0];
                        double *grad_y_so_z = gradient_so_z[3 * iatom + 1];
                        double *grad_z_so_z = gradient_so_z[3 * iatom + 2];

                        add_block_to_matrix(dim, dim, grad_x_so_z, bra_dim, ket_dim, buf_grad_so_z[0], ioffset, joffset,
                                            1.0);
                        add_block_to_matrix(dim, dim, grad_y_so_z, bra_dim, ket_dim, buf_grad_so_z[1], ioffset, joffset,
                                            1.0);
                        add_block_to_matrix(dim, dim, grad_z_so_z, bra_dim, ket_dim, buf_grad_so_z[2], ioffset, joffset,
                                            1.0);

                        dealloc_gradients(buf_grad_so_x);
                        dealloc_gradients(buf_grad_so_y);
                        dealloc_gradients(buf_grad_so_z);
                    }
                }
            }

            double t2 = abs_time();
            printf("%10.3f sec\n", t2 - t1);

            joffset += ket_dim;
        }

        ioffset += bra_dim;
    }
}


/**
 * This function evaluates gradients of grpp integrals with respect to coordinates
 * of all atoms in a given molecule.
 */
void evaluate_outercore_potential_integrals_gradient(
        int num_shells,
        libgrpp_shell_t **shell_list,
        molecule_t *molecule,
        libgrpp_grpp_t **grpp_list,
        double **gradient_arep,
        double **gradient_so_x,
        double **gradient_so_y,
        double **gradient_so_z
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
            printf("outercore pot grad: ishell=%3d (L=%d)\tjshell=%3d (L=%d)\t", ishell + 1, bra->L, jshell + 1, ket->L);

            for (int iatom = 0; iatom < molecule->n_atoms; iatom++) {

                double point_3d[3];
                point_3d[0] = molecule->coord_x[iatom];
                point_3d[1] = molecule->coord_y[iatom];
                point_3d[2] = molecule->coord_z[iatom];

                for (int irpp = 0; irpp < molecule->n_atoms; irpp++) {
                    int z = molecule->charges[irpp];
                    libgrpp_grpp_t *grpp_operator = grpp_list[z];
                    if (grpp_operator == NULL) {
                        continue;
                    }

                    double grpp_origin[3];
                    grpp_origin[0] = molecule->coord_x[irpp];
                    grpp_origin[1] = molecule->coord_y[irpp];
                    grpp_origin[2] = molecule->coord_z[irpp];

                    double **buf_grad_arep = alloc_gradients(bra, ket);
                    double **buf_grad_so_x = alloc_gradients(bra, ket);
                    double **buf_grad_so_y = alloc_gradients(bra, ket);
                    double **buf_grad_so_z = alloc_gradients(bra, ket);

                    libgrpp_outercore_potential_integrals_gradient(
                            bra, ket, grpp_origin,
                            grpp_operator->n_oc_shells, grpp_operator->U_oc, grpp_operator->oc_shells, point_3d,
                            buf_grad_arep, buf_grad_so_x, buf_grad_so_y, buf_grad_so_z
                            );

                    /*
                     * AREP
                     */
                    double *grad_x_arep = gradient_arep[3 * iatom + 0];
                    double *grad_y_arep = gradient_arep[3 * iatom + 1];
                    double *grad_z_arep = gradient_arep[3 * iatom + 2];

                    add_block_to_matrix(dim, dim, grad_x_arep, bra_dim, ket_dim, buf_grad_arep[0], ioffset, joffset,
                                        1.0);
                    add_block_to_matrix(dim, dim, grad_y_arep, bra_dim, ket_dim, buf_grad_arep[1], ioffset, joffset,
                                        1.0);
                    add_block_to_matrix(dim, dim, grad_z_arep, bra_dim, ket_dim, buf_grad_arep[2], ioffset, joffset,
                                        1.0);

                    /*
                     * SO-X
                     */
                    double *grad_x_so_x = gradient_so_x[3 * iatom + 0];
                    double *grad_y_so_x = gradient_so_x[3 * iatom + 1];
                    double *grad_z_so_x = gradient_so_x[3 * iatom + 2];

                    add_block_to_matrix(dim, dim, grad_x_so_x, bra_dim, ket_dim, buf_grad_so_x[0], ioffset, joffset,
                                        1.0);
                    add_block_to_matrix(dim, dim, grad_y_so_x, bra_dim, ket_dim, buf_grad_so_x[1], ioffset, joffset,
                                        1.0);
                    add_block_to_matrix(dim, dim, grad_z_so_x, bra_dim, ket_dim, buf_grad_so_x[2], ioffset, joffset,
                                        1.0);

                    /*
                     * SO-Y
                     */
                    double *grad_x_so_y = gradient_so_y[3 * iatom + 0];
                    double *grad_y_so_y = gradient_so_y[3 * iatom + 1];
                    double *grad_z_so_y = gradient_so_y[3 * iatom + 2];

                    add_block_to_matrix(dim, dim, grad_x_so_y, bra_dim, ket_dim, buf_grad_so_y[0], ioffset, joffset,
                                        1.0);
                    add_block_to_matrix(dim, dim, grad_y_so_y, bra_dim, ket_dim, buf_grad_so_y[1], ioffset, joffset,
                                        1.0);
                    add_block_to_matrix(dim, dim, grad_z_so_y, bra_dim, ket_dim, buf_grad_so_y[2], ioffset, joffset,
                                        1.0);

                    /*
                     * SO-Z
                     */
                    double *grad_x_so_z = gradient_so_z[3 * iatom + 0];
                    double *grad_y_so_z = gradient_so_z[3 * iatom + 1];
                    double *grad_z_so_z = gradient_so_z[3 * iatom + 2];

                    add_block_to_matrix(dim, dim, grad_x_so_z, bra_dim, ket_dim, buf_grad_so_z[0], ioffset, joffset,
                                        1.0);
                    add_block_to_matrix(dim, dim, grad_y_so_z, bra_dim, ket_dim, buf_grad_so_z[1], ioffset, joffset,
                                        1.0);
                    add_block_to_matrix(dim, dim, grad_z_so_z, bra_dim, ket_dim, buf_grad_so_z[2], ioffset, joffset,
                                        1.0);

                    dealloc_gradients(buf_grad_arep);
                    dealloc_gradients(buf_grad_so_x);
                    dealloc_gradients(buf_grad_so_y);
                    dealloc_gradients(buf_grad_so_z);
                }
            }

            double t2 = abs_time();
            printf("%10.3f sec\n", t2 - t1);

            joffset += ket_dim;
        }

        ioffset += bra_dim;
    }
}


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

