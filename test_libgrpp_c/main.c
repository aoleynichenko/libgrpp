/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2022 Alexander Oleynichenko
 */

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "../libgrpp/libgrpp.h"

#include "xyz.h"
#include "basis.h"
#include "rpp.h"
#include "shell_list.h"
#include "eval_integrals.h"
#include "overlap_gradients.h"
#include "grpp_gradients.h"
#include "print_matrix.h"
#include "abs_time.h"


int main(int argc, char **argv)
{
    basis_set_t *basis_sets[N_CHEM_ELEMENTS];
    grpp_t *grecps[N_CHEM_ELEMENTS];

    libgrpp_set_cartesian_order(LIBGRPP_CART_ORDER_TURBOMOLE);

    printf("\n");
    printf("    -----------------------------------------------------\n");
    printf("           the C front-end to the libgrecp library       \n");
    printf("    -----------------------------------------------------\n");
    printf("    a. oleynichenko                            3 dec 2021\n");
    printf("    -----------------------------------------------------\n");
    printf("\n");

    /* read molecular geometry */
    molecule_t *molecule = read_molecule("molecule.xyz");
    if (molecule == NULL) {
        return 1;
    }
    print_molecule(stdout, molecule);

    /* find basis sets for each type of atoms */
    memset(basis_sets, 0, sizeof(basis_sets));
    for (int iatom = 0; iatom < molecule->n_atoms; iatom++) {
        int nuc_charge = molecule->charges[iatom];
        if (basis_sets[nuc_charge] == NULL) {
            basis_set_t *set = read_basis_set("basis.inp", nuc_charge);
            basis_sets[nuc_charge] = set;
        }
    }

    /* print basis sets */
    for (int i = 0; i < N_CHEM_ELEMENTS; i++) {
        if (basis_sets[i] == NULL) {
            continue;
        }
        printf("\n\tbasis set for element Z = %d:\n", i);
        printf("\t------------------------------\n\n");
        print_basis_set(stdout, basis_sets[i]);
    }

    /* read ECPs for each type of atoms */
    memset(grecps, 0, sizeof(grecps));
    for (int iatom = 0; iatom < molecule->n_atoms; iatom++) {
        int nuc_charge = molecule->charges[iatom];
        if (grecps[nuc_charge] == NULL) {
            grpp_t *grecp = read_grpp("ecp.inp", nuc_charge);
            grecps[nuc_charge] = grecp;
        }
    }

    /* print ECPs */
    for (int i = 0; i < N_CHEM_ELEMENTS; i++) {
        if (grecps[i] == NULL) {
            continue;
        }
        printf("\n\teffective core potential for element Z = %d:\n", i);
        printf("\t---------------------------------------------\n\n");
        print_grpp(stdout, grecps[i]);
    }

    /* generate list of atom-centered basis shells */
    int num_shells;
    libgrpp_shell_t **shell_list = construct_shell_list(molecule, basis_sets, &num_shells);
    print_shell_list(shell_list, num_shells);

    /* ECP integral evaluation */
    int basis_dim = calculate_basis_dim(shell_list, num_shells);
    double *arep_matrix = (double *) calloc(basis_dim * basis_dim, sizeof(double));
    double *so_x_matrix = (double *) calloc(basis_dim * basis_dim, sizeof(double));
    double *so_y_matrix = (double *) calloc(basis_dim * basis_dim, sizeof(double));
    double *so_z_matrix = (double *) calloc(basis_dim * basis_dim, sizeof(double));
    double *nucattr_matrix = (double *) calloc(basis_dim * basis_dim, sizeof(double));
    double *overlap_matrix = (double *) calloc(basis_dim * basis_dim, sizeof(double));

    double time_start = abs_time();

    libgrpp_set_radial_tolerance(1e-16);
    evaluate_grpp_integrals(num_shells, shell_list, molecule, grecps,
                       arep_matrix, so_x_matrix, so_y_matrix, so_z_matrix);

    double time_finish = abs_time();
    printf("\nECP integration time: %.3f sec\n\n", time_finish - time_start);

    /*time_start = abs_time();
    evaluate_nuclear_attraction_integrals(num_shells, shell_list, molecule, nucattr_matrix, LIBGRPP_NUCLEAR_MODEL_FERMI);
    time_finish = abs_time();
    printf("\nNuclear attraction (Coulomb) integration time: %.3f sec\n\n", time_finish - time_start);*/

    //evaluate_overlap_integrals_gradients(num_shells, shell_list, molecule);
    //evaluate_grpp_integrals_gradients(num_shells, shell_list, molecule, grecps);

    /* print matrices to files */
    print_matrix_lower_triangle("libgrpp_c_arep.txt", basis_dim, arep_matrix);
    print_matrix_lower_triangle("libgrpp_c_so_x.txt", basis_dim, so_x_matrix);
    print_matrix_lower_triangle("libgrpp_c_so_y.txt", basis_dim, so_y_matrix);
    print_matrix_lower_triangle("libgrpp_c_so_z.txt", basis_dim, so_z_matrix);
    //print_matrix_lower_triangle("libgrpp_c_coulomb.txt", basis_dim, nucattr_matrix);
    //print_matrix_lower_triangle("libgrpp_c_overlap.txt", basis_dim, overlap_matrix);

    /* cleanup */
    free(arep_matrix);
    free(so_x_matrix);
    free(so_y_matrix);
    free(so_z_matrix);
    free(nucattr_matrix);
    free(overlap_matrix);

    delete_molecule(molecule);
    delete_shell_list(shell_list, num_shells);
    for (int i = 0; i < N_CHEM_ELEMENTS; i++) {
        if (basis_sets[i] != NULL) {
            delete_basis_set(basis_sets[i]);
        }
        if (grecps[i] != NULL) {
            delete_grpp(grecps[i]);
        }
    }
}
