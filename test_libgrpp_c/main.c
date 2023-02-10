/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2023 Alexander Oleynichenko
 */

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "../libgrpp/libgrpp.h"

#include "molecule.h"
#include "basis.h"
#include "rpp.h"
#include "shell_list.h"
#include "eval_integrals.h"
#include "eval_integrals_grad.h"
#include "grpp_gradients.h"
#include "print_matrix.h"
#include "abs_time.h"


int main(int argc, char **argv)
{
    basis_set_t *basis_sets[N_CHEM_ELEMENTS];
    libgrpp_grpp_t *grpps[N_CHEM_ELEMENTS];

    // turn off buffering for stdout (for seeing output immediately)
    setvbuf(stdout, NULL, _IONBF, 0);

    libgrpp_set_cartesian_order(LIBGRPP_CART_ORDER_TURBOMOLE);

    printf("\n");
    printf("    -----------------------------------------------------\n");
    printf("           the C front-end to the libgrpp library        \n");
    printf("    -----------------------------------------------------\n");
    printf("    a. oleynichenko                            9 feb 2023\n");
    printf("    -----------------------------------------------------\n");
    printf("\n");

    /*
     * read molecular geometry
     */
    molecule_t *molecule = read_molecule("molecule.xyz");
    if (molecule == NULL) {
        return 1;
    }
    print_molecule(stdout, molecule);

    /*
     * find basis sets for each type of atoms
     */
    memset(basis_sets, 0, sizeof(basis_sets));
    for (int iatom = 0; iatom < molecule->n_atoms; iatom++) {
        int nuc_charge = molecule->charges[iatom];
        if (basis_sets[nuc_charge] == NULL) {
            basis_set_t *set = read_basis_set("basis.inp", nuc_charge);
            basis_sets[nuc_charge] = set;
        }
    }

    /*
     * print basis sets
     */
    for (int i = 0; i < N_CHEM_ELEMENTS; i++) {
        if (basis_sets[i] == NULL) {
            continue;
        }
        printf("\n\tbasis set for element Z = %d:\n", i);
        printf("\t------------------------------\n\n");
        print_basis_set(stdout, basis_sets[i]);
    }

    /*
     * read ECPs for each type of atoms
     */
    memset(grpps, 0, sizeof(grpps));
    for (int iatom = 0; iatom < molecule->n_atoms; iatom++) {
        int nuc_charge = molecule->charges[iatom];
        if (grpps[nuc_charge] == NULL) {
            libgrpp_grpp_t *grecp = read_grpp("ecp.inp", nuc_charge);
            grpps[nuc_charge] = grecp;
        }
    }

    /*
     * print ECPs
     */
    for (int i = 0; i < N_CHEM_ELEMENTS; i++) {
        if (grpps[i] == NULL) {
            continue;
        }
        printf("\n\teffective core potential for element Z = %d:\n", i);
        printf("\t---------------------------------------------\n\n");
        print_grpp(stdout, grpps[i]);
    }

    /*
     * generate list of atom-centered basis shells
     */
    int num_shells;
    libgrpp_shell_t **shell_list = construct_shell_list(molecule, basis_sets, &num_shells);
    print_shell_list(shell_list, num_shells);

    /*
     * pseudopotential integrals
     */
    double time_start = abs_time();

    int basis_dim = calculate_basis_dim(shell_list, num_shells);
    double *arep_matrix = (double *) calloc(basis_dim * basis_dim, sizeof(double));
    double *so_x_matrix = (double *) calloc(basis_dim * basis_dim, sizeof(double));
    double *so_y_matrix = (double *) calloc(basis_dim * basis_dim, sizeof(double));
    double *so_z_matrix = (double *) calloc(basis_dim * basis_dim, sizeof(double));

    libgrpp_set_radial_tolerance(1e-16);
    evaluate_grpp_integrals(num_shells, shell_list, molecule, grpps,
                            arep_matrix, so_x_matrix, so_y_matrix, so_z_matrix);

    double time_finish = abs_time();
    printf("\ntime for pseudopotential integrals: %.3f sec\n\n", time_finish - time_start);

    /*
     * overlap integrals
     */
    time_start = abs_time();
    double *overlap_matrix = (double *) calloc(basis_dim * basis_dim, sizeof(double));

    evaluate_overlap_integrals(num_shells, shell_list, overlap_matrix);

    time_finish = abs_time();
    printf("\ntime for overlap integrals: %.3f sec\n\n", time_finish - time_start);

    /*
     * nuclear attraction integrals
     */
    time_start = abs_time();
    double *nucattr_matrix = (double *) calloc(basis_dim * basis_dim, sizeof(double));

    evaluate_nuclear_attraction_integrals(num_shells, shell_list, molecule, nucattr_matrix, LIBGRPP_NUCLEAR_MODEL_POINT_CHARGE);

    time_finish = abs_time();
    printf("\ntime for nuclear attraction integrals: %.3f sec\n\n", time_finish - time_start);

    /*
     * print matrices to files
     */
    print_matrix_lower_triangle("libgrpp_c_arep.txt", basis_dim, arep_matrix);
    print_matrix_lower_triangle("libgrpp_c_so_x.txt", basis_dim, so_x_matrix);
    print_matrix_lower_triangle("libgrpp_c_so_y.txt", basis_dim, so_y_matrix);
    print_matrix_lower_triangle("libgrpp_c_so_z.txt", basis_dim, so_z_matrix);
    print_matrix_lower_triangle("libgrpp_c_overlap.txt", basis_dim, overlap_matrix);
    print_matrix_lower_triangle("libgrpp_c_nucattr.txt", basis_dim, nucattr_matrix);

    /*
     * gradients with respect to nuclear coordinates: overlap integrals
     */
    time_start = abs_time();

    double **grad = (double **) calloc(3 * molecule->n_atoms, sizeof(double *));
    for (int icoord = 0; icoord < 3 * molecule->n_atoms; icoord++) {
        grad[icoord] = (double *) calloc(basis_dim * basis_dim, sizeof(double));
    }

    evaluate_overlap_integrals_gradient(num_shells, shell_list, molecule, grad);
    for (int iatom = 0; iatom < molecule->n_atoms; iatom++) {
        char file_name_buf_x[100];
        char file_name_buf_y[100];
        char file_name_buf_z[100];

        sprintf(file_name_buf_x, "libgrpp_c_overlap_grad_%dx.txt", iatom);
        sprintf(file_name_buf_y, "libgrpp_c_overlap_grad_%dy.txt", iatom);
        sprintf(file_name_buf_z, "libgrpp_c_overlap_grad_%dz.txt", iatom);

        print_matrix_lower_triangle(file_name_buf_x, basis_dim, grad[3 * iatom + 0]);
        print_matrix_lower_triangle(file_name_buf_y, basis_dim, grad[3 * iatom + 1]);
        print_matrix_lower_triangle(file_name_buf_z, basis_dim, grad[3 * iatom + 2]);
    }

    time_finish = abs_time();
    printf("\ntime for overlap integrals gradients: %.3f sec\n\n", time_finish - time_start);

    //evaluate_grpp_integrals_gradients(num_shells, shell_list, molecule, grpps);

    /*
     * cleanup
     */
    free(arep_matrix);
    free(so_x_matrix);
    free(so_y_matrix);
    free(so_z_matrix);
    free(nucattr_matrix);
    free(overlap_matrix);

    for (int icoord = 0; icoord < 3 * molecule->n_atoms; icoord++) {
        free(grad[icoord]);
    }
    free(grad);

    delete_molecule(molecule);
    delete_shell_list(shell_list, num_shells);
    for (int i = 0; i < N_CHEM_ELEMENTS; i++) {
        if (basis_sets[i] != NULL) {
            delete_basis_set(basis_sets[i]);
        }
        if (grpps[i] != NULL) {
            libgrpp_delete_grpp(grpps[i]);
        }
    }
}
