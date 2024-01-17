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
#include "print_matrix.h"
#include "abs_time.h"


void calculate_write_grpp_integrals(int num_shells, libgrpp_shell_t **shell_list, molecule_t *molecule,
                                    libgrpp_grpp_t **grpps);

void calculate_write_overlap_integrals(int num_shells, libgrpp_shell_t **shell_list);

void calculate_write_nuclear_attraction_integrals(int num_shells, libgrpp_shell_t **shell_list, molecule_t *molecule,
                                                  int nuclear_model);

void calculate_write_overlap_gradient(int num_shells, libgrpp_shell_t **shell_list, molecule_t *molecule);

void calculate_write_kinetic_energy_integrals(int num_shells, libgrpp_shell_t **shell_list);

void calculate_write_grpp_gradient(int num_shells, libgrpp_shell_t **shell_list, molecule_t *molecule,
                                   libgrpp_grpp_t **grpps);

void calculate_write_type1_gradient(int num_shells, libgrpp_shell_t **shell_list, molecule_t *molecule,
                                   libgrpp_grpp_t **grpps);

void calculate_write_type2_gradient(int num_shells, libgrpp_shell_t **shell_list, molecule_t *molecule,
                                    libgrpp_grpp_t **grpps);

void calculate_write_spin_orbit_gradient(int num_shells, libgrpp_shell_t **shell_list, molecule_t *molecule,
                                         libgrpp_grpp_t **grpps);

void calculate_write_outercore_potential_gradient(int num_shells, libgrpp_shell_t **shell_list, molecule_t *molecule,
                                                  libgrpp_grpp_t **grpps);

void calculate_write_momentum_integrals(int num_shells, libgrpp_shell_t **shell_list);


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
    printf("    a. oleynichenko                           25 dec 2023\n");
    printf("    -----------------------------------------------------\n");
    printf("\n");

    /*
     * parse command-line arguments
     */
    int calc_grpp_integrals = 0;
    int calc_ovlp_integrals = 0;
    int calc_coul_point_integrals = 0;
    int calc_coul_point_num_integrals = 0;
    int calc_coul_ball_integrals = 0;
    int calc_coul_gauss_integrals = 0;
    int calc_coul_fermi_integrals = 0;
    int calc_coul_fermi_bubble_integrals = 0;
    int calc_kine_integrals = 0;
    int calc_mome_integrals = 0;
    int calc_ovlp_gradients = 0;
    int calc_grpp_gradients = 0;
    int calc_type1_gradients = 0;
    int calc_type2_gradients = 0;
    int calc_spin_orbit_gradients = 0;
    int calc_outercore_gradients = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--grpp") == 0) {
            calc_grpp_integrals = 1;
        }
        else if (strcmp(argv[i], "--overlap") == 0) {
            calc_ovlp_integrals = 1;
        }
        else if (strcmp(argv[i], "--coulomb-point") == 0) {
            calc_coul_point_integrals = 1;
        }
        else if (strcmp(argv[i], "--coulomb-point-num") == 0) {
            calc_coul_point_num_integrals = 1;
        }
        else if (strcmp(argv[i], "--coulomb-ball") == 0) {
            calc_coul_ball_integrals = 1;
        }
        else if (strcmp(argv[i], "--coulomb-gauss") == 0) {
            calc_coul_gauss_integrals = 1;
        }
        else if (strcmp(argv[i], "--coulomb-fermi") == 0) {
            calc_coul_fermi_integrals = 1;
        }
        else if (strcmp(argv[i], "--coulomb-fermi-bubble") == 0) {
            calc_coul_fermi_bubble_integrals = 1;
        }
        else if (strcmp(argv[i], "--kinetic") == 0) {
            calc_kine_integrals = 1;
        }
        else if (strcmp(argv[i], "--momentum") == 0) {
            calc_mome_integrals = 1;
        }
        else if (strcmp(argv[i], "--grpp-grad") == 0) {
            calc_grpp_gradients = 1;
        }
        else if (strcmp(argv[i], "--grad-type1") == 0) {
            calc_type1_gradients = 1;
        }
        else if (strcmp(argv[i], "--grad-type2") == 0) {
            calc_type2_gradients = 1;
        }
        else if (strcmp(argv[i], "--grad-spin-orbit") == 0) {
            calc_spin_orbit_gradients = 1;
        }
        else if (strcmp(argv[i], "--grad-outercore") == 0) {
            calc_outercore_gradients = 1;
        }
        else if (strcmp(argv[i], "--overlap-grad") == 0) {
            calc_ovlp_gradients = 1;
        }
        else {
            printf("unknown argument: %s\n", argv[1]);
            printf("\n");
            printf("list of possible arguments (types of integrals to be calculated):\n");
            printf("--grpp                   generalized pseudopotential\n");
            printf("--overlap                overlap\n");
            printf("--coulomb-point          nuclear attraction, point nucleus\n");
            printf("--coulomb-point-num      nuclear attraction, point nucleus (numerical radial integrals)\n");
            printf("--coulomb-ball           nuclear attraction, charged ball distribution\n");
            printf("--coulomb-gauss          nuclear attraction, gaussian distribution\n");
            printf("--coulomb-fermi          nuclear attraction, fermi distribution\n");
            printf("--coulomb-fermi-bubble   nuclear attraction, fermi distribution with bubble\n");
            printf("--kinetic                kinetic energy operator\n");
            printf("--momentum               momentum operator\n");
            printf("--grpp-grad              gradients of generalized pseudopotential integrals\n");
            printf("--grad-type1             gradients of the type1 pseudopotential integrals\n");
            printf("--grad-type2             gradients of the type2 pseudopotential integrals\n");
            printf("--grad-spin-orbit        gradients of spin-orbit (pseudopotential) integrals\n");
            printf("--grad-outercore         gradients of outercore potential (pseudopotential) integrals\n");
            printf("--overlap-grad           gradients of overlap integrals\n");
            return 1;
        }
    }

    printf("\n");
    printf(" types of integrals to be calculated:\n");
    printf("\n");
    printf(" pseudopotential (grpp)      %s\n", calc_grpp_integrals ? "yes" : "no");
    printf(" coulomb (point nuc)         %s\n", calc_coul_point_integrals ? "yes" : "no");
    printf(" coulomb (point nuc, num)    %s\n", calc_coul_point_num_integrals ? "yes" : "no");
    printf(" coulomb (charged ball nuc)  %s\n", calc_coul_ball_integrals ? "yes" : "no");
    printf(" coulomb (gauss nuc)         %s\n", calc_coul_gauss_integrals ? "yes" : "no");
    printf(" coulomb (fermi nuc)         %s\n", calc_coul_fermi_integrals ? "yes" : "no");
    printf(" coulomb (fermi bubble nuc)  %s\n", calc_coul_fermi_bubble_integrals ? "yes" : "no");
    printf(" overlap                     %s\n", calc_ovlp_integrals ? "yes" : "no");
    printf(" kinetic-energy              %s\n", calc_kine_integrals ? "yes" : "no");
    printf(" momentum                    %s\n", calc_mome_integrals ? "yes" : "no");
    printf(" grpp gradients              %s\n", calc_grpp_gradients ? "yes" : "no");
    printf(" pp type1 gradients          %s\n", calc_type1_gradients ? "yes" : "no");
    printf(" pp type1 gradients          %s\n", calc_type2_gradients ? "yes" : "no");
    printf(" pp spin-orbit gradients     %s\n", calc_spin_orbit_gradients ? "yes" : "no");
    printf(" pp outercore pot gradients  %s\n", calc_outercore_gradients ? "yes" : "no");
    printf(" overlap gradients           %s\n", calc_ovlp_gradients ? "yes" : "no");
    printf("\n");

    /*
     * initialize LIBGRPP
     */
    libgrpp_init();

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
    if (calc_grpp_integrals) {
        calculate_write_grpp_integrals(num_shells, shell_list, molecule, grpps);
    }

    /*
     * overlap integrals
     */
    if (calc_ovlp_integrals) {
        calculate_write_overlap_integrals(num_shells, shell_list);
    }

    /*
     * kinetic energy integrals
     */
    if (calc_kine_integrals) {
        calculate_write_kinetic_energy_integrals(num_shells, shell_list);
    }


    /*
     * momentum integrals (imaginary part)
     */
    if (calc_mome_integrals) {
        calculate_write_momentum_integrals(num_shells, shell_list);
    }


    /*
     * nuclear attraction integrals
     * (for different models of nuclear charge distribution)
     */
    if (calc_coul_point_integrals) {
        calculate_write_nuclear_attraction_integrals(num_shells, shell_list, molecule,
                                                     LIBGRPP_NUCLEAR_MODEL_POINT_CHARGE);
    }

    if (calc_coul_point_num_integrals) {
        calculate_write_nuclear_attraction_integrals(num_shells, shell_list, molecule,
                                                     LIBGRPP_NUCLEAR_MODEL_POINT_CHARGE_NUMERICAL);
    }

    if (calc_coul_ball_integrals) {
        calculate_write_nuclear_attraction_integrals(num_shells, shell_list, molecule,
                                                     LIBGRPP_NUCLEAR_MODEL_CHARGED_BALL);
    }

    if (calc_coul_gauss_integrals) {
        calculate_write_nuclear_attraction_integrals(num_shells, shell_list, molecule, LIBGRPP_NUCLEAR_MODEL_GAUSSIAN);
    }

    if (calc_coul_fermi_integrals) {
        calculate_write_nuclear_attraction_integrals(num_shells, shell_list, molecule, LIBGRPP_NUCLEAR_MODEL_FERMI);
    }

    if (calc_coul_fermi_bubble_integrals) {
        calculate_write_nuclear_attraction_integrals(num_shells, shell_list, molecule,
                                                     LIBGRPP_NUCLEAR_MODEL_FERMI_BUBBLE);
    }

    /*
     * gradients with respect to nuclear coordinates: overlap integrals
     */
    if (calc_ovlp_gradients) {
        calculate_write_overlap_gradient(num_shells, shell_list, molecule);
    }

    /*
     * gradients with respect to nuclear coordinates: GRPP integrals
     */
    if (calc_grpp_gradients) {
        calculate_write_grpp_gradient(num_shells, shell_list, molecule, grpps);
    }

    if (calc_type1_gradients) {
        calculate_write_type1_gradient(num_shells, shell_list, molecule, grpps);
    }

    if (calc_type2_gradients) {
        calculate_write_type2_gradient(num_shells, shell_list, molecule, grpps);
    }

    if (calc_spin_orbit_gradients) {
        calculate_write_spin_orbit_gradient(num_shells, shell_list, molecule, grpps);
    }

    if (calc_outercore_gradients) {
        calculate_write_outercore_potential_gradient(num_shells, shell_list, molecule, grpps);
    }

    /*
     * cleanup
     */
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


void calculate_write_grpp_integrals(int num_shells, libgrpp_shell_t **shell_list, molecule_t *molecule,
                                    libgrpp_grpp_t **grpps)
{
    int basis_dim = calculate_basis_dim(shell_list, num_shells);

    double *arep_matrix = (double *) calloc(basis_dim * basis_dim, sizeof(double));
    double *so_x_matrix = (double *) calloc(basis_dim * basis_dim, sizeof(double));
    double *so_y_matrix = (double *) calloc(basis_dim * basis_dim, sizeof(double));
    double *so_z_matrix = (double *) calloc(basis_dim * basis_dim, sizeof(double));

    double time_start = abs_time();

    libgrpp_set_radial_tolerance(1e-16);
    evaluate_grpp_integrals(num_shells, shell_list, molecule, grpps,
                            arep_matrix, so_x_matrix, so_y_matrix, so_z_matrix);

    double time_finish = abs_time();

    print_matrix_lower_triangle("libgrpp_c_arep.txt", basis_dim, arep_matrix);
    print_matrix_lower_triangle("libgrpp_c_so_x.txt", basis_dim, so_x_matrix);
    print_matrix_lower_triangle("libgrpp_c_so_y.txt", basis_dim, so_y_matrix);
    print_matrix_lower_triangle("libgrpp_c_so_z.txt", basis_dim, so_z_matrix);

    free(arep_matrix);
    free(so_x_matrix);
    free(so_y_matrix);
    free(so_z_matrix);

    printf("\ntime for pseudopotential integrals: %.3f sec\n\n", time_finish - time_start);
}


void calculate_write_overlap_integrals(int num_shells, libgrpp_shell_t **shell_list)
{
    int basis_dim = calculate_basis_dim(shell_list, num_shells);
    double *overlap_matrix = (double *) calloc(basis_dim * basis_dim, sizeof(double));

    double time_start = abs_time();
    evaluate_overlap_integrals(num_shells, shell_list, overlap_matrix);
    double time_finish = abs_time();

    print_matrix_lower_triangle("libgrpp_c_overlap.txt", basis_dim, overlap_matrix);

    free(overlap_matrix);

    printf("\ntime for overlap integrals: %.3f sec\n\n", time_finish - time_start);
}


void calculate_write_kinetic_energy_integrals(int num_shells, libgrpp_shell_t **shell_list)
{
    int basis_dim = calculate_basis_dim(shell_list, num_shells);
    double *kinetic_matrix = (double *) calloc(basis_dim * basis_dim, sizeof(double));

    double time_start = abs_time();
    evaluate_kinetic_energy_integrals(num_shells, shell_list, kinetic_matrix);
    double time_finish = abs_time();

    print_matrix_lower_triangle("libgrpp_c_kinetic_energy.txt", basis_dim, kinetic_matrix);

    free(kinetic_matrix);

    printf("\ntime for kinetic energy integrals: %.3f sec\n\n", time_finish - time_start);
}


void calculate_write_momentum_integrals(int num_shells, libgrpp_shell_t **shell_list)
{
    int basis_dim = calculate_basis_dim(shell_list, num_shells);
    double *px_matrix = (double *) calloc(basis_dim * basis_dim, sizeof(double));
    double *py_matrix = (double *) calloc(basis_dim * basis_dim, sizeof(double));
    double *pz_matrix = (double *) calloc(basis_dim * basis_dim, sizeof(double));

    double time_start = abs_time();
    evaluate_momentum_integrals(num_shells, shell_list, px_matrix, py_matrix, pz_matrix);
    double time_finish = abs_time();

    print_matrix_lower_triangle("libgrpp_c_momentum_x.txt", basis_dim, px_matrix);
    print_matrix_lower_triangle("libgrpp_c_momentum_y.txt", basis_dim, py_matrix);
    print_matrix_lower_triangle("libgrpp_c_momentum_z.txt", basis_dim, pz_matrix);

    free(px_matrix);
    free(py_matrix);
    free(pz_matrix);

    printf("\ntime for momentum integrals: %.3f sec\n\n", time_finish - time_start);
}


void calculate_write_nuclear_attraction_integrals(int num_shells, libgrpp_shell_t **shell_list, molecule_t *molecule,
                                                  int nuclear_model)
{
    int basis_dim = calculate_basis_dim(shell_list, num_shells);
    double *nucattr_matrix = (double *) calloc(basis_dim * basis_dim, sizeof(double));

    double time_start = abs_time();
    evaluate_nuclear_attraction_integrals(num_shells, shell_list, molecule, nucattr_matrix, nuclear_model);
    double time_finish = abs_time();

    switch (nuclear_model) {
        case LIBGRPP_NUCLEAR_MODEL_POINT_CHARGE:
            print_matrix_lower_triangle("libgrpp_c_nucattr_point.txt", basis_dim, nucattr_matrix);
            break;
        case LIBGRPP_NUCLEAR_MODEL_POINT_CHARGE_NUMERICAL:
            print_matrix_lower_triangle("libgrpp_c_nucattr_point_numerical.txt", basis_dim, nucattr_matrix);
            break;
        case LIBGRPP_NUCLEAR_MODEL_CHARGED_BALL:
            print_matrix_lower_triangle("libgrpp_c_nucattr_charged_ball.txt", basis_dim, nucattr_matrix);
            break;
        case LIBGRPP_NUCLEAR_MODEL_GAUSSIAN:
            print_matrix_lower_triangle("libgrpp_c_nucattr_gaussian.txt", basis_dim, nucattr_matrix);
            break;
        case LIBGRPP_NUCLEAR_MODEL_FERMI:
            print_matrix_lower_triangle("libgrpp_c_nucattr_fermi.txt", basis_dim, nucattr_matrix);
            break;
        case LIBGRPP_NUCLEAR_MODEL_FERMI_BUBBLE:
            print_matrix_lower_triangle("libgrpp_c_nucattr_fermi_bubble.txt", basis_dim, nucattr_matrix);
            break;
        default:
            break;
    }

    free(nucattr_matrix);

    printf("\ntime for nuclear attraction integrals: %.3f sec\n\n", time_finish - time_start);
}


void calculate_write_overlap_gradient(int num_shells, libgrpp_shell_t **shell_list, molecule_t *molecule)
{
    int basis_dim = calculate_basis_dim(shell_list, num_shells);

    double **grad = (double **) calloc(3 * molecule->n_atoms, sizeof(double *));
    for (int icoord = 0; icoord < 3 * molecule->n_atoms; icoord++) {
        grad[icoord] = (double *) calloc(basis_dim * basis_dim, sizeof(double));
    }

    double time_start = abs_time();
    evaluate_overlap_integrals_gradient(num_shells, shell_list, molecule, grad);
    double time_finish = abs_time();

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

    for (int icoord = 0; icoord < 3 * molecule->n_atoms; icoord++) {
        free(grad[icoord]);
    }
    free(grad);

    printf("\ntime for overlap integrals gradients: %.3f sec\n\n", time_finish - time_start);
}


void calculate_write_grpp_gradient(int num_shells, libgrpp_shell_t **shell_list, molecule_t *molecule,
                                   libgrpp_grpp_t **grpps)
{
    int basis_dim = calculate_basis_dim(shell_list, num_shells);

    double **grad_arep = (double **) calloc(3 * molecule->n_atoms, sizeof(double *));
    double **grad_so_x = (double **) calloc(3 * molecule->n_atoms, sizeof(double *));
    double **grad_so_y = (double **) calloc(3 * molecule->n_atoms, sizeof(double *));
    double **grad_so_z = (double **) calloc(3 * molecule->n_atoms, sizeof(double *));

    for (int icoord = 0; icoord < 3 * molecule->n_atoms; icoord++) {
        grad_arep[icoord] = (double *) calloc(basis_dim * basis_dim, sizeof(double));
        grad_so_x[icoord] = (double *) calloc(basis_dim * basis_dim, sizeof(double));
        grad_so_y[icoord] = (double *) calloc(basis_dim * basis_dim, sizeof(double));
        grad_so_z[icoord] = (double *) calloc(basis_dim * basis_dim, sizeof(double));
    }

    double time_start = abs_time();
    evaluate_grpp_integrals_gradient(num_shells, shell_list, molecule, grpps, grad_arep, grad_so_x, grad_so_y,
                                     grad_so_z);
    double time_finish = abs_time();
    printf("\ntime for grpp integrals gradients: %.3f sec\n\n", time_finish - time_start);

    for (int iatom = 0; iatom < molecule->n_atoms; iatom++) {

        /*
         * AREP gradients
         */
        char file_name_buf_arep_x[100];
        char file_name_buf_arep_y[100];
        char file_name_buf_arep_z[100];

        sprintf(file_name_buf_arep_x, "libgrpp_c_arep_grad_%dx.txt", iatom);
        sprintf(file_name_buf_arep_y, "libgrpp_c_arep_grad_%dy.txt", iatom);
        sprintf(file_name_buf_arep_z, "libgrpp_c_arep_grad_%dz.txt", iatom);

        print_matrix_lower_triangle(file_name_buf_arep_x, basis_dim, grad_arep[3 * iatom + 0]);
        print_matrix_lower_triangle(file_name_buf_arep_y, basis_dim, grad_arep[3 * iatom + 1]);
        print_matrix_lower_triangle(file_name_buf_arep_z, basis_dim, grad_arep[3 * iatom + 2]);

        /*
         * SO-X gradients
         */
        char file_name_buf_so_x_x[100];
        char file_name_buf_so_x_y[100];
        char file_name_buf_so_x_z[100];

        sprintf(file_name_buf_so_x_x, "libgrpp_c_so_x_grad_%dx.txt", iatom);
        sprintf(file_name_buf_so_x_y, "libgrpp_c_so_x_grad_%dy.txt", iatom);
        sprintf(file_name_buf_so_x_z, "libgrpp_c_so_x_grad_%dz.txt", iatom);

        print_matrix_lower_triangle(file_name_buf_so_x_x, basis_dim, grad_so_x[3 * iatom + 0]);
        print_matrix_lower_triangle(file_name_buf_so_x_y, basis_dim, grad_so_x[3 * iatom + 1]);
        print_matrix_lower_triangle(file_name_buf_so_x_z, basis_dim, grad_so_x[3 * iatom + 2]);

        /*
         * SO-Y gradients
         */
        char file_name_buf_so_y_x[100];
        char file_name_buf_so_y_y[100];
        char file_name_buf_so_y_z[100];

        sprintf(file_name_buf_so_y_x, "libgrpp_c_so_y_grad_%dx.txt", iatom);
        sprintf(file_name_buf_so_y_y, "libgrpp_c_so_y_grad_%dy.txt", iatom);
        sprintf(file_name_buf_so_y_z, "libgrpp_c_so_y_grad_%dz.txt", iatom);

        print_matrix_lower_triangle(file_name_buf_so_y_x, basis_dim, grad_so_y[3 * iatom + 0]);
        print_matrix_lower_triangle(file_name_buf_so_y_y, basis_dim, grad_so_y[3 * iatom + 1]);
        print_matrix_lower_triangle(file_name_buf_so_y_z, basis_dim, grad_so_y[3 * iatom + 2]);

        /*
         * SO-Z gradients
         */
        char file_name_buf_so_z_x[100];
        char file_name_buf_so_z_y[100];
        char file_name_buf_so_z_z[100];

        sprintf(file_name_buf_so_z_x, "libgrpp_c_so_z_grad_%dx.txt", iatom);
        sprintf(file_name_buf_so_z_y, "libgrpp_c_so_z_grad_%dy.txt", iatom);
        sprintf(file_name_buf_so_z_z, "libgrpp_c_so_z_grad_%dz.txt", iatom);

        print_matrix_lower_triangle(file_name_buf_so_z_x, basis_dim, grad_so_z[3 * iatom + 0]);
        print_matrix_lower_triangle(file_name_buf_so_z_y, basis_dim, grad_so_z[3 * iatom + 1]);
        print_matrix_lower_triangle(file_name_buf_so_z_z, basis_dim, grad_so_z[3 * iatom + 2]);
    }

    /*
     * cleanup
     */
    for (int icoord = 0; icoord < 3 * molecule->n_atoms; icoord++) {
        free(grad_arep[icoord]);
        free(grad_so_x[icoord]);
        free(grad_so_y[icoord]);
        free(grad_so_z[icoord]);
    }
    free(grad_arep);
    free(grad_so_x);
    free(grad_so_y);
    free(grad_so_z);
}


void calculate_write_type1_gradient(int num_shells, libgrpp_shell_t **shell_list, molecule_t *molecule,
                                    libgrpp_grpp_t **grpps)
{
    int basis_dim = calculate_basis_dim(shell_list, num_shells);

    double **grad_arep = (double **) calloc(3 * molecule->n_atoms, sizeof(double *));
    for (int icoord = 0; icoord < 3 * molecule->n_atoms; icoord++) {
        grad_arep[icoord] = (double *) calloc(basis_dim * basis_dim, sizeof(double));
    }

    double time_start = abs_time();
    evaluate_type1_integrals_gradient(num_shells, shell_list, molecule, grpps, grad_arep);
    double time_finish = abs_time();
    printf("\ntime for type1 pp integrals gradients: %.3f sec\n\n", time_finish - time_start);

    for (int iatom = 0; iatom < molecule->n_atoms; iatom++) {
        /*
         * AREP gradients
         */
        char file_name_buf_arep_x[100];
        char file_name_buf_arep_y[100];
        char file_name_buf_arep_z[100];

        sprintf(file_name_buf_arep_x, "libgrpp_c_arep_type1_grad_%dx.txt", iatom);
        sprintf(file_name_buf_arep_y, "libgrpp_c_arep_type1_grad_%dy.txt", iatom);
        sprintf(file_name_buf_arep_z, "libgrpp_c_arep_type1_grad_%dz.txt", iatom);

        print_matrix_lower_triangle(file_name_buf_arep_x, basis_dim, grad_arep[3 * iatom + 0]);
        print_matrix_lower_triangle(file_name_buf_arep_y, basis_dim, grad_arep[3 * iatom + 1]);
        print_matrix_lower_triangle(file_name_buf_arep_z, basis_dim, grad_arep[3 * iatom + 2]);
    }

    /*
     * cleanup
     */
    for (int icoord = 0; icoord < 3 * molecule->n_atoms; icoord++) {
        free(grad_arep[icoord]);
    }
    free(grad_arep);
}


void calculate_write_type2_gradient(int num_shells, libgrpp_shell_t **shell_list, molecule_t *molecule,
                                    libgrpp_grpp_t **grpps)
{
    int basis_dim = calculate_basis_dim(shell_list, num_shells);

    double **grad_arep = (double **) calloc(3 * molecule->n_atoms, sizeof(double *));
    for (int icoord = 0; icoord < 3 * molecule->n_atoms; icoord++) {
        grad_arep[icoord] = (double *) calloc(basis_dim * basis_dim, sizeof(double));
    }

    double time_start = abs_time();
    evaluate_type2_integrals_gradient(num_shells, shell_list, molecule, grpps, grad_arep);
    double time_finish = abs_time();
    printf("\ntime for type2 pp integrals gradients: %.3f sec\n\n", time_finish - time_start);

    for (int iatom = 0; iatom < molecule->n_atoms; iatom++) {
        /*
         * AREP gradients
         */
        char file_name_buf_arep_x[100];
        char file_name_buf_arep_y[100];
        char file_name_buf_arep_z[100];

        sprintf(file_name_buf_arep_x, "libgrpp_c_arep_type2_grad_%dx.txt", iatom);
        sprintf(file_name_buf_arep_y, "libgrpp_c_arep_type2_grad_%dy.txt", iatom);
        sprintf(file_name_buf_arep_z, "libgrpp_c_arep_type2_grad_%dz.txt", iatom);

        print_matrix_lower_triangle(file_name_buf_arep_x, basis_dim, grad_arep[3 * iatom + 0]);
        print_matrix_lower_triangle(file_name_buf_arep_y, basis_dim, grad_arep[3 * iatom + 1]);
        print_matrix_lower_triangle(file_name_buf_arep_z, basis_dim, grad_arep[3 * iatom + 2]);
    }

    /*
     * cleanup
     */
    for (int icoord = 0; icoord < 3 * molecule->n_atoms; icoord++) {
        free(grad_arep[icoord]);
    }
    free(grad_arep);
}


void calculate_write_spin_orbit_gradient(int num_shells, libgrpp_shell_t **shell_list, molecule_t *molecule,
                                   libgrpp_grpp_t **grpps)
{
    int basis_dim = calculate_basis_dim(shell_list, num_shells);

    double **grad_so_x = (double **) calloc(3 * molecule->n_atoms, sizeof(double *));
    double **grad_so_y = (double **) calloc(3 * molecule->n_atoms, sizeof(double *));
    double **grad_so_z = (double **) calloc(3 * molecule->n_atoms, sizeof(double *));

    for (int icoord = 0; icoord < 3 * molecule->n_atoms; icoord++) {
        grad_so_x[icoord] = (double *) calloc(basis_dim * basis_dim, sizeof(double));
        grad_so_y[icoord] = (double *) calloc(basis_dim * basis_dim, sizeof(double));
        grad_so_z[icoord] = (double *) calloc(basis_dim * basis_dim, sizeof(double));
    }

    double time_start = abs_time();
    evaluate_spin_orbit_integrals_gradient(num_shells, shell_list, molecule, grpps, grad_so_x, grad_so_y,
                                     grad_so_z);
    double time_finish = abs_time();
    printf("\ntime for spin-orbit integrals gradients: %.3f sec\n\n", time_finish - time_start);

    for (int iatom = 0; iatom < molecule->n_atoms; iatom++) {

        /*
         * SO-X gradients
         */
        char file_name_buf_so_x_x[100];
        char file_name_buf_so_x_y[100];
        char file_name_buf_so_x_z[100];

        sprintf(file_name_buf_so_x_x, "libgrpp_c_so_x_so_grad_%dx.txt", iatom);
        sprintf(file_name_buf_so_x_y, "libgrpp_c_so_x_so_grad_%dy.txt", iatom);
        sprintf(file_name_buf_so_x_z, "libgrpp_c_so_x_so_grad_%dz.txt", iatom);

        print_matrix_lower_triangle(file_name_buf_so_x_x, basis_dim, grad_so_x[3 * iatom + 0]);
        print_matrix_lower_triangle(file_name_buf_so_x_y, basis_dim, grad_so_x[3 * iatom + 1]);
        print_matrix_lower_triangle(file_name_buf_so_x_z, basis_dim, grad_so_x[3 * iatom + 2]);

        /*
         * SO-Y gradients
         */
        char file_name_buf_so_y_x[100];
        char file_name_buf_so_y_y[100];
        char file_name_buf_so_y_z[100];

        sprintf(file_name_buf_so_y_x, "libgrpp_c_so_y_so_grad_%dx.txt", iatom);
        sprintf(file_name_buf_so_y_y, "libgrpp_c_so_y_so_grad_%dy.txt", iatom);
        sprintf(file_name_buf_so_y_z, "libgrpp_c_so_y_so_grad_%dz.txt", iatom);

        print_matrix_lower_triangle(file_name_buf_so_y_x, basis_dim, grad_so_y[3 * iatom + 0]);
        print_matrix_lower_triangle(file_name_buf_so_y_y, basis_dim, grad_so_y[3 * iatom + 1]);
        print_matrix_lower_triangle(file_name_buf_so_y_z, basis_dim, grad_so_y[3 * iatom + 2]);

        /*
         * SO-Z gradients
         */
        char file_name_buf_so_z_x[100];
        char file_name_buf_so_z_y[100];
        char file_name_buf_so_z_z[100];

        sprintf(file_name_buf_so_z_x, "libgrpp_c_so_z_so_grad_%dx.txt", iatom);
        sprintf(file_name_buf_so_z_y, "libgrpp_c_so_z_so_grad_%dy.txt", iatom);
        sprintf(file_name_buf_so_z_z, "libgrpp_c_so_z_so_grad_%dz.txt", iatom);

        print_matrix_lower_triangle(file_name_buf_so_z_x, basis_dim, grad_so_z[3 * iatom + 0]);
        print_matrix_lower_triangle(file_name_buf_so_z_y, basis_dim, grad_so_z[3 * iatom + 1]);
        print_matrix_lower_triangle(file_name_buf_so_z_z, basis_dim, grad_so_z[3 * iatom + 2]);
    }

    /*
     * cleanup
     */
    for (int icoord = 0; icoord < 3 * molecule->n_atoms; icoord++) {
        free(grad_so_x[icoord]);
        free(grad_so_y[icoord]);
        free(grad_so_z[icoord]);
    }
    free(grad_so_x);
    free(grad_so_y);
    free(grad_so_z);
}


void calculate_write_outercore_potential_gradient(int num_shells, libgrpp_shell_t **shell_list, molecule_t *molecule,
                                   libgrpp_grpp_t **grpps)
{
    int basis_dim = calculate_basis_dim(shell_list, num_shells);

    double **grad_arep = (double **) calloc(3 * molecule->n_atoms, sizeof(double *));
    double **grad_so_x = (double **) calloc(3 * molecule->n_atoms, sizeof(double *));
    double **grad_so_y = (double **) calloc(3 * molecule->n_atoms, sizeof(double *));
    double **grad_so_z = (double **) calloc(3 * molecule->n_atoms, sizeof(double *));

    for (int icoord = 0; icoord < 3 * molecule->n_atoms; icoord++) {
        grad_arep[icoord] = (double *) calloc(basis_dim * basis_dim, sizeof(double));
        grad_so_x[icoord] = (double *) calloc(basis_dim * basis_dim, sizeof(double));
        grad_so_y[icoord] = (double *) calloc(basis_dim * basis_dim, sizeof(double));
        grad_so_z[icoord] = (double *) calloc(basis_dim * basis_dim, sizeof(double));
    }

    double time_start = abs_time();
    evaluate_outercore_potential_integrals_gradient(
            num_shells, shell_list, molecule, grpps,
            grad_arep, grad_so_x, grad_so_y, grad_so_z
    );

    double time_finish = abs_time();
    printf("\ntime for outercore potential integrals gradients: %.3f sec\n\n", time_finish - time_start);

    for (int iatom = 0; iatom < molecule->n_atoms; iatom++) {

        /*
         * AREP gradients
         */
        char file_name_buf_arep_x[100];
        char file_name_buf_arep_y[100];
        char file_name_buf_arep_z[100];

        sprintf(file_name_buf_arep_x, "libgrpp_c_arep_ocpot_grad_%dx.txt", iatom);
        sprintf(file_name_buf_arep_y, "libgrpp_c_arep_ocpot_grad_%dy.txt", iatom);
        sprintf(file_name_buf_arep_z, "libgrpp_c_arep_ocpot_grad_%dz.txt", iatom);

        print_matrix_lower_triangle(file_name_buf_arep_x, basis_dim, grad_arep[3 * iatom + 0]);
        print_matrix_lower_triangle(file_name_buf_arep_y, basis_dim, grad_arep[3 * iatom + 1]);
        print_matrix_lower_triangle(file_name_buf_arep_z, basis_dim, grad_arep[3 * iatom + 2]);

        /*
         * SO-X gradients
         */
        char file_name_buf_so_x_x[100];
        char file_name_buf_so_x_y[100];
        char file_name_buf_so_x_z[100];

        sprintf(file_name_buf_so_x_x, "libgrpp_c_so_x_ocpot_grad_%dx.txt", iatom);
        sprintf(file_name_buf_so_x_y, "libgrpp_c_so_x_ocpot_grad_%dy.txt", iatom);
        sprintf(file_name_buf_so_x_z, "libgrpp_c_so_x_ocpot_grad_%dz.txt", iatom);

        print_matrix_lower_triangle(file_name_buf_so_x_x, basis_dim, grad_so_x[3 * iatom + 0]);
        print_matrix_lower_triangle(file_name_buf_so_x_y, basis_dim, grad_so_x[3 * iatom + 1]);
        print_matrix_lower_triangle(file_name_buf_so_x_z, basis_dim, grad_so_x[3 * iatom + 2]);

        /*
         * SO-Y gradients
         */
        char file_name_buf_so_y_x[100];
        char file_name_buf_so_y_y[100];
        char file_name_buf_so_y_z[100];

        sprintf(file_name_buf_so_y_x, "libgrpp_c_so_y_ocpot_grad_%dx.txt", iatom);
        sprintf(file_name_buf_so_y_y, "libgrpp_c_so_y_ocpot_grad_%dy.txt", iatom);
        sprintf(file_name_buf_so_y_z, "libgrpp_c_so_y_ocpot_grad_%dz.txt", iatom);

        print_matrix_lower_triangle(file_name_buf_so_y_x, basis_dim, grad_so_y[3 * iatom + 0]);
        print_matrix_lower_triangle(file_name_buf_so_y_y, basis_dim, grad_so_y[3 * iatom + 1]);
        print_matrix_lower_triangle(file_name_buf_so_y_z, basis_dim, grad_so_y[3 * iatom + 2]);

        /*
         * SO-Z gradients
         */
        char file_name_buf_so_z_x[100];
        char file_name_buf_so_z_y[100];
        char file_name_buf_so_z_z[100];

        sprintf(file_name_buf_so_z_x, "libgrpp_c_so_z_ocpot_grad_%dx.txt", iatom);
        sprintf(file_name_buf_so_z_y, "libgrpp_c_so_z_ocpot_grad_%dy.txt", iatom);
        sprintf(file_name_buf_so_z_z, "libgrpp_c_so_z_ocpot_grad_%dz.txt", iatom);

        print_matrix_lower_triangle(file_name_buf_so_z_x, basis_dim, grad_so_z[3 * iatom + 0]);
        print_matrix_lower_triangle(file_name_buf_so_z_y, basis_dim, grad_so_z[3 * iatom + 1]);
        print_matrix_lower_triangle(file_name_buf_so_z_z, basis_dim, grad_so_z[3 * iatom + 2]);
    }

    /*
     * cleanup
     */
    for (int icoord = 0; icoord < 3 * molecule->n_atoms; icoord++) {
        free(grad_arep[icoord]);
        free(grad_so_x[icoord]);
        free(grad_so_y[icoord]);
        free(grad_so_z[icoord]);
    }
    free(grad_arep);
    free(grad_so_x);
    free(grad_so_y);
    free(grad_so_z);
}
