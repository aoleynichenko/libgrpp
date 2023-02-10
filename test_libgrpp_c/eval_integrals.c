/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2023 Alexander Oleynichenko
 */

#include "eval_integrals.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "abs_time.h"
#include "shell_list.h"

#define MAX_BUF 10000


void update_vector(int size, double *a, double factor, double *b);

void add_block_to_matrix(int dim_1, int dim_2, double *matrix,
                         int block_dim_1, int block_dim_2, double *block,
                         int col_offset, int row_offset, double factor);

void evaluate_grpp_integrals_shell_pair(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        libgrpp_grpp_t *grpp_operator,
        double *grpp_origin,
        double *arep_matrix,
        double *so_x_matrix,
        double *so_y_matrix,
        double *so_z_matrix
);


/**
 * evaluates matrix elements of the GRPP operator
 */
void evaluate_grpp_integrals(int num_shells, libgrpp_shell_t **shell_list,
                             molecule_t *molecule, libgrpp_grpp_t **grpp_list,
                             double *arep_matrix, double *so_x_matrix, double *so_y_matrix, double *so_z_matrix)
{
    double buf_arep[MAX_BUF];
    double buf_spin_orbit[3][MAX_BUF];

    int dim = calculate_basis_dim(shell_list, num_shells);

    memset(arep_matrix, 0, sizeof(double) * dim * dim);
    memset(so_x_matrix, 0, sizeof(double) * dim * dim);
    memset(so_y_matrix, 0, sizeof(double) * dim * dim);
    memset(so_z_matrix, 0, sizeof(double) * dim * dim);

    int ioffset = 0;
    for (int ishell = 0; ishell < num_shells; ishell++) {

        libgrpp_shell_t *bra = shell_list[ishell];
        int bra_dim = libgrpp_get_shell_size(bra);

        int joffset = 0;
        for (int jshell = 0; jshell < num_shells; jshell++) {

            libgrpp_shell_t *ket = shell_list[jshell];
            int ket_dim = libgrpp_get_shell_size(ket);

            double t1 = abs_time();
            printf("grpp: ishell=%3d (L=%d)\tjshell=%3d (L=%d)\t", ishell + 1, bra->L, jshell + 1, ket->L);

            for (int iatom = 0; iatom < molecule->n_atoms; iatom++) {
                int z = molecule->charges[iatom];
                libgrpp_grpp_t *grpp = grpp_list[z];
                if (grpp == NULL) {
                    continue;
                }

                double ecp_origin[3];
                ecp_origin[0] = molecule->coord_x[iatom];
                ecp_origin[1] = molecule->coord_y[iatom];
                ecp_origin[2] = molecule->coord_z[iatom];

                evaluate_grpp_integrals_shell_pair(
                        bra, ket, grpp, ecp_origin,
                        buf_arep, buf_spin_orbit[0], buf_spin_orbit[1], buf_spin_orbit[2]
                );

                add_block_to_matrix(dim, dim, arep_matrix, bra_dim, ket_dim, buf_arep, ioffset, joffset, 1.0);
                add_block_to_matrix(dim, dim, so_x_matrix, bra_dim, ket_dim, buf_spin_orbit[0], ioffset, joffset, 1.0);
                add_block_to_matrix(dim, dim, so_y_matrix, bra_dim, ket_dim, buf_spin_orbit[1], ioffset, joffset, 1.0);
                add_block_to_matrix(dim, dim, so_z_matrix, bra_dim, ket_dim, buf_spin_orbit[2], ioffset, joffset, 1.0);
            }

            double t2 = abs_time();
            printf("%10.3f sec\n", t2 - t1);

            joffset += ket_dim;
        }

        ioffset += bra_dim;
    }
}


void evaluate_grpp_integrals_shell_pair(
        libgrpp_shell_t *shell_A,
        libgrpp_shell_t *shell_B,
        libgrpp_grpp_t *grpp_operator,
        double *grpp_origin,
        double *arep_matrix,
        double *so_x_matrix,
        double *so_y_matrix,
        double *so_z_matrix
)
{
    size_t size = shell_A->cart_size * shell_B->cart_size;
    double *buf_arep = (double *) calloc(size, sizeof(double));
    double *buf_so_x = (double *) calloc(size, sizeof(double));
    double *buf_so_y = (double *) calloc(size, sizeof(double));
    double *buf_so_z = (double *) calloc(size, sizeof(double));

    memset(arep_matrix, 0, sizeof(double) * size);
    memset(so_x_matrix, 0, sizeof(double) * size);
    memset(so_y_matrix, 0, sizeof(double) * size);
    memset(so_z_matrix, 0, sizeof(double) * size);

    /*
     * radially-local ("type-1") integrals
     */
    libgrpp_type1_integrals(shell_A, shell_B, grpp_origin, grpp_operator->U_L, buf_arep);
    update_vector(size, arep_matrix, 1.0, buf_arep);

    /*
     * semilocal AREP ("type-2") integrals
     */
    for (int L = 0; L < grpp_operator->n_arep; L++) {
        libgrpp_type2_integrals(shell_A, shell_B, grpp_origin, grpp_operator->U_arep[L], buf_arep);
        update_vector(size, arep_matrix, 1.0, buf_arep);
    }

    /*
     * semilocal SO ("type-3") integrals
     */
    for (int L = 1; L < grpp_operator->n_esop; L++) {
        libgrpp_spin_orbit_integrals(shell_A, shell_B, grpp_origin, grpp_operator->U_esop[L],
                                     buf_so_x, buf_so_y, buf_so_z);

        update_vector(size, so_x_matrix, 2.0 / (2 * L + 1), buf_so_x);
        update_vector(size, so_y_matrix, 2.0 / (2 * L + 1), buf_so_y);
        update_vector(size, so_z_matrix, 2.0 / (2 * L + 1), buf_so_z);
    }

    /*
     * integrals over outercore non-local potentials,
     * the part specific for GRPP.
     *
     * note that proper pre-factors for the SO part are calculated inside
     * the libgrpp_outercore_potential_integrals() procedure.
     */
    libgrpp_outercore_potential_integrals(
            shell_A, shell_B, grpp_origin,
            grpp_operator->n_oc_shells, grpp_operator->U_oc, grpp_operator->oc_shells,
            buf_arep, buf_so_x, buf_so_y, buf_so_z
    );

    update_vector(size, arep_matrix, 1.0, buf_arep);
    update_vector(size, so_x_matrix, 1.0, buf_so_x);
    update_vector(size, so_y_matrix, 1.0, buf_so_y);
    update_vector(size, so_z_matrix, 1.0, buf_so_z);

    /*
     * cleanup
     */
    free(buf_arep);
    free(buf_so_x);
    free(buf_so_y);
    free(buf_so_z);
}


/**
 * evaluates overlap matrix
 */
void evaluate_overlap_integrals(int num_shells, libgrpp_shell_t **shell_list, double *overlap_matrix)
{
    double buf[MAX_BUF];

    int dim = calculate_basis_dim(shell_list, num_shells);

    memset(overlap_matrix, 0, sizeof(double) * dim * dim);

    int ioffset = 0;
    for (int ishell = 0; ishell < num_shells; ishell++) {

        libgrpp_shell_t *bra_shell = shell_list[ishell];
        int bra_dim = libgrpp_get_shell_size(bra_shell);

        int joffset = 0;
        for (int jshell = 0; jshell < num_shells; jshell++) {

            libgrpp_shell_t *ket_shell = shell_list[jshell];
            int ket_dim = libgrpp_get_shell_size(ket_shell);

            double t1 = abs_time();
            printf("overlap: ishell=%3d (L=%d)\tjshell=%3d (L=%d)\t", ishell + 1, bra_shell->L, jshell + 1, ket_shell->L);

            libgrpp_overlap_integrals(bra_shell, ket_shell, buf);
            add_block_to_matrix(dim, dim, overlap_matrix, bra_dim, ket_dim, buf, ioffset, joffset, 1.0);

            double t2 = abs_time();
            printf("%10.3f sec\n", t2 - t1);

            joffset += ket_dim;
        }

        ioffset += bra_dim;
    }
}


/**
 * evaluates nuclear attraction integrals
 */
void evaluate_nuclear_attraction_integrals(int num_shells, libgrpp_shell_t **shell_list, molecule_t *molecule,
                                           double *nucattr_matrix, int nuclear_model)
{
    double buf_nucattr[MAX_BUF];

    int dim = calculate_basis_dim(shell_list, num_shells);

    memset(nucattr_matrix, 0, sizeof(double) * dim * dim);

    int ioffset = 0;
    for (int ishell = 0; ishell < num_shells; ishell++) {

        libgrpp_shell_t *bra = shell_list[ishell];
        int bra_dim = libgrpp_get_shell_size(bra);

        int joffset = 0;
        for (int jshell = 0; jshell < num_shells; jshell++) {

            libgrpp_shell_t *ket = shell_list[jshell];
            int ket_dim = libgrpp_get_shell_size(ket);

            double t1 = abs_time();
            printf("nuc attr: ishell=%3d (L=%d)\tjshell=%3d (L=%d)\t", ishell + 1, bra->L, jshell + 1, ket->L);

            for (int iatom = 0; iatom < molecule->n_atoms; iatom++) {
                double charge_origin[3];
                charge_origin[0] = molecule->coord_x[iatom];
                charge_origin[1] = molecule->coord_y[iatom];
                charge_origin[2] = molecule->coord_z[iatom];
                int z = molecule->charges[iatom];

                double nucmod_params[10];
                int mass_number = get_element_mass_number_abundant(z);
                double R_rms = libgrpp_estimate_nuclear_rms_radius_johnson_1985(mass_number);

                if (nuclear_model == LIBGRPP_NUCLEAR_MODEL_POINT_CHARGE) {
                    // do nothing
                }
                else if (nuclear_model == LIBGRPP_NUCLEAR_MODEL_GAUSSIAN ||
                         nuclear_model == LIBGRPP_NUCLEAR_MODEL_CHARGED_BALL) {
                    nucmod_params[0] = R_rms * FERMI_UNITS_TO_ATOMIC;
                }
                else { // Fermi
                    double c, a;
                    libgrpp_estimate_fermi_model_parameters(R_rms, &c, &a);
                    nucmod_params[0] = c * FERMI_UNITS_TO_ATOMIC;
                    nucmod_params[1] = a * FERMI_UNITS_TO_ATOMIC;
                }

                libgrpp_nuclear_attraction_integrals(bra, ket, charge_origin, z, nuclear_model, nucmod_params,
                                                     buf_nucattr);

                add_block_to_matrix(dim, dim, nucattr_matrix, bra_dim, ket_dim, buf_nucattr, ioffset, joffset, 1.0);
            }

            double t2 = abs_time();
            printf("%10.3f sec\n", t2 - t1);

            joffset += ket_dim;
        }

        ioffset += bra_dim;
    }
}


void add_block_to_matrix(int dim_1, int dim_2, double *matrix,
                         int block_dim_1, int block_dim_2, double *block,
                         int col_offset, int row_offset, double factor)
{
    for (int i = 0; i < block_dim_1; i++) {
        for (int j = 0; j < block_dim_2; j++) {
            double element = block[i * block_dim_2 + j];
            matrix[(col_offset + i) * dim_2 + (row_offset + j)] += factor * element;
        }
    }
}


/**
 * a += factor * b
 */
void update_vector(int size, double *a, double factor, double *b)
{
    for (int i = 0; i < size; i++) {
        a[i] += factor * b[i];
    }
}

