/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2023 Alexander Oleynichenko
 */

/*
 * Wrappers for the LIBGRPP subroutines to be used from Fortran projects.
 *
 * C99       Fortran
 * --------------------
 * int32_t   integer(4)
 * double    real(8)
 */

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "factorial.h"
#include "libgrpp.h"

/*
 * Fine-tuning of the LIBGRPP internal parameters.
 */

void libgrpp_set_default_parameters_()
{
    libgrpp_set_default_parameters();
}


void libgrpp_set_radial_tolerance_(const double *tolerance)
{
    libgrpp_set_radial_tolerance(*tolerance);
}


void libgrpp_set_angular_screening_tolerance_(const double *tolerance)
{
    libgrpp_set_angular_screening_tolerance(*tolerance);
}


void libgrpp_set_modified_bessel_tolerance_(const double *tolerance)
{
    libgrpp_set_modified_bessel_tolerance(*tolerance);
}


void libgrpp_set_cartesian_order_(const int32_t *order)
{
    libgrpp_set_cartesian_order(*order);
}


/*
 * Type 1 RPP integrals (local term)
 */

void libgrpp_type1_integrals_(
        // contracted Gaussian A
        double *origin_A,
        int32_t *L_A,
        int32_t *num_primitives_A,
        double *coeffs_A,
        double *alpha_A,
        // contracted Gaussian B
        double *origin_B,
        int32_t *L_B,
        int32_t *num_primitives_B,
        double *coeffs_B,
        double *alpha_B,
        // pseudopotential
        double *rpp_origin,
        int32_t *rpp_num_primitives,
        int32_t *rpp_powers,
        double *rpp_coeffs,
        double *rpp_alpha,
        // answer
        double *matrix
)
{
    int *pot_powers_int = (int *) calloc(*rpp_num_primitives, sizeof(int));

    for (int i = 0; i < *rpp_num_primitives; i++) {
        pot_powers_int[i] = rpp_powers[i];
    }

    libgrpp_potential_t *pot = libgrpp_new_potential(0, 0, *rpp_num_primitives, pot_powers_int, rpp_coeffs, rpp_alpha);
    libgrpp_shell_t *shell_A = libgrpp_new_shell(origin_A, *L_A, *num_primitives_A, coeffs_A, alpha_A);
    libgrpp_shell_t *shell_B = libgrpp_new_shell(origin_B, *L_B, *num_primitives_B, coeffs_B, alpha_B);

    libgrpp_type1_integrals(shell_A, shell_B, rpp_origin, pot, matrix);

    libgrpp_delete_shell(shell_A);
    libgrpp_delete_shell(shell_B);
    libgrpp_delete_potential(pot);
    free(pot_powers_int);
}


/*
 * Type 2 RPP integrals (semilocal terms with projectors)
 */

void libgrpp_type2_integrals_(
        // contracted Gaussian A
        double *origin_A,
        int32_t *L_A,
        int32_t *num_primitives_A,
        double *coeffs_A,
        double *alpha_A,
        double *origin_B,
        // contracted Gaussian B
        int32_t *L_B,
        int32_t *num_primitives_B,
        double *coeffs_B,
        double *alpha_B,
        // pseudopotential
        double *pot_origin,
        int32_t *pot_L,
        int32_t *pot_num_primitives,
        int32_t *pot_powers,
        double *pot_coeffs,
        double *pot_alpha,
        // answer
        double *matrix
)
{
    int *pot_powers_int = (int *) calloc(*pot_num_primitives, sizeof(int));

    for (int i = 0; i < *pot_num_primitives; i++) {
        pot_powers_int[i] = pot_powers[i];
    }

    libgrpp_potential_t *pot = libgrpp_new_potential(*pot_L, 0, *pot_num_primitives, pot_powers_int, pot_coeffs,
                                                     pot_alpha);
    libgrpp_shell_t *shell_A = libgrpp_new_shell(origin_A, *L_A, *num_primitives_A, coeffs_A, alpha_A);
    libgrpp_shell_t *shell_B = libgrpp_new_shell(origin_B, *L_B, *num_primitives_B, coeffs_B, alpha_B);

    libgrpp_type2_integrals(shell_A, shell_B, pot_origin, pot, matrix);

    libgrpp_delete_shell(shell_A);
    libgrpp_delete_shell(shell_B);
    libgrpp_delete_potential(pot);
    free(pot_powers_int);
}


/*
 * Effective spin-orbit operator ("Type 3") RPP integrals (semilocal terms with projectors)
 */

void libgrpp_spin_orbit_integrals_(
        // contracted Gaussian A
        double *origin_A,
        int32_t *L_A,
        int32_t *num_primitives_A,
        double *coeffs_A,
        double *alpha_A,
        // contracted Gaussian B
        double *origin_B,
        int32_t *L_B,
        int32_t *num_primitives_B,
        double *coeffs_B,
        double *alpha_B,
        // pseudopotential
        double *pot_origin,
        int32_t *pot_L,
        int32_t *pot_num_primitives,
        int32_t *pot_powers,
        double *pot_coeffs,
        double *pot_alpha,
        // answer
        double *so_x_matrix,
        double *so_y_matrix,
        double *so_z_matrix
)
{
    int *pot_powers_int = (int *) calloc(*pot_num_primitives, sizeof(int));

    for (int i = 0; i < *pot_num_primitives; i++) {
        pot_powers_int[i] = pot_powers[i];
    }

    /*
     * construct RPP structure
     */
    libgrpp_potential_t *pot = libgrpp_new_potential(*pot_L, 0, *pot_num_primitives, pot_powers_int, pot_coeffs,
                                                     pot_alpha);

    /*
     * construct shells
     */
    libgrpp_shell_t *shell_A = libgrpp_new_shell(origin_A, *L_A, *num_primitives_A, coeffs_A, alpha_A);
    libgrpp_shell_t *shell_B = libgrpp_new_shell(origin_B, *L_B, *num_primitives_B, coeffs_B, alpha_B);

    libgrpp_spin_orbit_integrals(shell_A, shell_B, pot_origin, pot, so_x_matrix, so_y_matrix, so_z_matrix);

    libgrpp_delete_shell(shell_A);
    libgrpp_delete_shell(shell_B);
    libgrpp_delete_potential(pot);
    free(pot_powers_int);
}


/**
 * Outercore RPP integrals (non-local terms with projectors onto outercore spinors)
 *
 * Part 1: integration of the first non-local term:
 * U*|nlj><nlj| + |nlj><nlj|*U
 */
void libgrpp_outercore_potential_integrals_part_1_(
        // contracted Gaussian A
        double *origin_A,
        int32_t *L_A,
        int32_t *num_primitives_A,
        double *coeffs_A,
        double *alpha_A,
        double *origin_B,
        // contracted Gaussian B
        int32_t *L_B,
        int32_t *num_primitives_B,
        double *coeffs_B,
        double *alpha_B,
        // pseudopotential for the outercore shell LJ
        double *pot_origin,
        int32_t *pot_L,
        int32_t *pot_J,
        int32_t *pot_num_primitives,
        int32_t *pot_powers,
        double *pot_coeffs,
        double *pot_alpha,
        // expansion of the outercore shell LJ
        int32_t *oc_shell_num_primitives,
        double *oc_shell_coeffs,
        double *oc_shell_alpha,
        // answer
        double *arep_matrix,
        double *so_x_matrix,
        double *so_y_matrix,
        double *so_z_matrix
)
{
    void libgrpp_outercore_potential_integrals_part_1(
            libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B,
            double *C, libgrpp_potential_t *oc_potential, libgrpp_shell_t *oc_shell,
            double *arep_matrix, double *so_x_matrix, double *so_y_matrix, double *so_z_matrix
    );

    /*
     * array conversion: Fortran -> C
     */
    int *pot_powers_int = (int *) calloc(*pot_num_primitives, sizeof(int));
    for (int i = 0; i < *pot_num_primitives; i++) {
        pot_powers_int[i] = pot_powers[i];
    }

    /*
     * construct shells
     */
    libgrpp_shell_t *shell_A = libgrpp_new_shell(origin_A, *L_A, *num_primitives_A, coeffs_A, alpha_A);
    libgrpp_shell_t *shell_B = libgrpp_new_shell(origin_B, *L_B, *num_primitives_B, coeffs_B, alpha_B);

    /*
     * construct pseudopotential for the given L,J numbers
     */
    libgrpp_potential_t *oc_potential = libgrpp_new_potential(
            *pot_L, *pot_J, *pot_num_primitives,
            pot_powers_int, pot_coeffs, pot_alpha
    );

    /*
     * construct outercore shell associated with the pseudopotential
     */
    libgrpp_shell_t *oc_shell = libgrpp_new_shell(
            pot_origin, *pot_L, *oc_shell_num_primitives,
            oc_shell_coeffs, oc_shell_alpha
    );

    /*
     * evaluate RPP integrals
     */
    libgrpp_outercore_potential_integrals_part_1(
            shell_A, shell_B,
            pot_origin, oc_potential, oc_shell,
            arep_matrix, so_x_matrix, so_y_matrix, so_z_matrix
    );

    /*
     * clean-up
     */
    libgrpp_delete_shell(shell_A);
    libgrpp_delete_shell(shell_B);
    libgrpp_delete_potential(oc_potential);
    libgrpp_delete_shell(oc_shell);
    free(pot_powers_int);
}


/**
 * Outercore RPP integrals (non-local terms with projectors onto outercore spinors)
 *
 * Part 2: integration of the second non-local term:
 * |nlj><nlj| U |n'lj><n'lj|
 */
void libgrpp_outercore_potential_integrals_part_2_(
        // contracted Gaussian A
        double *origin_A,
        int32_t *L_A,
        int32_t *num_primitives_A,
        double *coeffs_A,
        double *alpha_A,
        // contracted Gaussian B
        double *origin_B,
        int32_t *L_B,
        int32_t *num_primitives_B,
        double *coeffs_B,
        double *alpha_B,
        // origin of the RPP
        double *pot_origin,
        // outercore shell 1:
        int32_t *oc_shell_1_L,
        int32_t *oc_shell_1_J,
        int32_t *pot1_num_primitives,
        int32_t *pot1_powers,
        double *pot1_coeffs,
        double *pot1_alpha,
        int32_t *oc_shell_1_num_primitives,
        double *oc_shell_1_coeffs,
        double *oc_shell_1_alpha,
        // outercore shell 2:
        int32_t *oc_shell_2_L,
        int32_t *oc_shell_2_J,
        int32_t *pot2_num_primitives,
        int32_t *pot2_powers,
        double *pot2_coeffs,
        double *pot2_alpha,
        int32_t *oc_shell_2_num_primitives,
        double *oc_shell_2_coeffs,
        double *oc_shell_2_alpha,
        // answer:
        double *arep_matrix,
        double *so_x_matrix,
        double *so_y_matrix,
        double *so_z_matrix
)
{
    void libgrpp_outercore_potential_integrals_part_2(
            libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B,
            double *C, libgrpp_potential_t *oc_potential_1, libgrpp_shell_t *oc_shell_1,
            libgrpp_potential_t *oc_potential_2, libgrpp_shell_t *oc_shell_2,
            double *arep_matrix, double *so_x_matrix, double *so_y_matrix, double *so_z_matrix
    );

    /*
     * array conversion: Fortran -> C
     */
    int *pot1_powers_int = (int *) calloc(*pot1_num_primitives, sizeof(int));
    int *pot2_powers_int = (int *) calloc(*pot2_num_primitives, sizeof(int));

    for (int i = 0; i < *pot1_num_primitives; i++) {
        pot1_powers_int[i] = pot1_powers[i];
    }
    for (int i = 0; i < *pot2_num_primitives; i++) {
        pot2_powers_int[i] = pot2_powers[i];
    }

    /*
     * construct shells
     */
    libgrpp_shell_t *shell_A = libgrpp_new_shell(origin_A, *L_A, *num_primitives_A, coeffs_A, alpha_A);
    libgrpp_shell_t *shell_B = libgrpp_new_shell(origin_B, *L_B, *num_primitives_B, coeffs_B, alpha_B);

    /*
     * the first outercore pseudopotential U_{n,L,J} and the corresponding outercore shell
     */
    libgrpp_potential_t *oc_potential_1 = libgrpp_new_potential(
            *oc_shell_1_L, *oc_shell_1_J, *pot1_num_primitives,
            pot1_powers_int, pot1_coeffs, pot1_alpha
    );
    libgrpp_shell_t *oc_shell_1 = libgrpp_new_shell(
            pot_origin, *oc_shell_1_L, *oc_shell_1_num_primitives,
            oc_shell_1_coeffs, oc_shell_1_alpha
    );

    /*
     * the second outercore pseudopotential U_{n',L',J'} and the corresponding outercore shell
     */
    libgrpp_potential_t *oc_potential_2 = libgrpp_new_potential(
            *oc_shell_2_L, *oc_shell_2_J, *pot2_num_primitives,
            pot2_powers_int, pot2_coeffs, pot2_alpha
    );
    libgrpp_shell_t *oc_shell_2 = libgrpp_new_shell(
            pot_origin, *oc_shell_2_L, *oc_shell_2_num_primitives,
            oc_shell_2_coeffs, oc_shell_2_alpha
    );

    /*
     * evaluate integrals
     */
    libgrpp_outercore_potential_integrals_part_2(
            shell_A, shell_B, pot_origin, oc_potential_1, oc_shell_1, oc_potential_2, oc_shell_2,
            arep_matrix, so_x_matrix, so_y_matrix, so_z_matrix
    );

    /*
     * clean-up
     */
    libgrpp_delete_shell(shell_A);
    libgrpp_delete_shell(shell_B);
    libgrpp_delete_potential(oc_potential_1);
    libgrpp_delete_potential(oc_potential_2);
    libgrpp_delete_shell(oc_shell_1);
    libgrpp_delete_shell(oc_shell_2);
    free(pot1_powers_int);
    free(pot2_powers_int);
}


/**
 * Overlap integrals between two contracted Gaussians with given cartesian parts x^n y^l z^m
 */

void evaluate_overlap_integral_contracted_(
        // contracted Gaussian A
        double *origin_A,
        int32_t *n_A,
        int32_t *l_A,
        int32_t *m_A,
        int32_t *num_primitives_A,
        double *coeffs_A,
        double *alpha_A,
        // contracted Gaussian B
        double *origin_B,
        int32_t *n_B,
        int32_t *l_B,
        int32_t *m_B,
        int32_t *num_primitives_B,
        double *coeffs_B,
        double *alpha_B,
        // answer
        double *overlap_integral
)
{
    void
    libgrpp_overlap_integrals(libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B, double *overlap_matrix);

    libgrpp_shell_t *shell_A = libgrpp_new_shell(origin_A, *n_A + *l_A + *m_A, *num_primitives_A, coeffs_A, alpha_A);
    libgrpp_shell_t *shell_B = libgrpp_new_shell(origin_B, *n_B + *l_B + *m_B, *num_primitives_B, coeffs_B, alpha_B);

    shell_A->cart_size = 1;
    shell_A->cart_list[0] = *n_A;
    shell_A->cart_list[1] = *l_A;
    shell_A->cart_list[2] = *m_A;

    shell_B->cart_size = 1;
    shell_B->cart_list[0] = *n_B;
    shell_B->cart_list[1] = *l_B;
    shell_B->cart_list[2] = *m_B;

    libgrpp_overlap_integrals(shell_A, shell_B, overlap_integral);

    libgrpp_delete_shell(shell_A);
    libgrpp_delete_shell(shell_B);
}


/*
 * calculates normalization factor for the given contracted Gaussians
 */
void radial_gto_norm_factor_(int32_t *L, int32_t *num_primitives, double *coeffs, double *alpha, double *norm)
{
    *norm = 0.0;
    double S = 0.0;
    double origin[] = {0, 0, 0};
    int zero = 0;

    evaluate_overlap_integral_contracted_(
            origin, L, &zero, &zero, num_primitives, coeffs, alpha,
            origin, L, &zero, &zero, num_primitives, coeffs, alpha,
            &S);

    *norm = sqrt(double_factorial(2 * (*L) - 1)) / sqrt(S);
}


/*
 * overlap integrals (for the shell pair)
 */

void libgrpp_overlap_integrals_(
        // contracted Gaussian A
        double *origin_A,
        int32_t *L_A,
        int32_t *num_primitives_A,
        double *coeffs_A,
        double *alpha_A,
        // contracted Gaussian B
        double *origin_B,
        int32_t *L_B,
        int32_t *num_primitives_B,
        double *coeffs_B,
        double *alpha_B,
        // answer
        double *matrix
)
{
    libgrpp_shell_t *shell_A = libgrpp_new_shell(origin_A, *L_A, *num_primitives_A, coeffs_A, alpha_A);
    libgrpp_shell_t *shell_B = libgrpp_new_shell(origin_B, *L_B, *num_primitives_B, coeffs_B, alpha_B);

    libgrpp_overlap_integrals(shell_A, shell_B, matrix);

    libgrpp_delete_shell(shell_A);
    libgrpp_delete_shell(shell_B);
}


/*
 * nuclear attraction integrals
 */

void libgrpp_nuclear_attraction_integrals_(
        // contracted Gaussian A
        double *origin_A,
        int32_t *L_A,
        int32_t *num_primitives_A,
        double *coeffs_A,
        double *alpha_A,
        // contracted Gaussian B
        double *origin_B,
        int32_t *L_B,
        int32_t *num_primitives_B,
        double *coeffs_B,
        double *alpha_B,
        // potential definition
        double *charge_origin,
        int32_t *charge,
        int32_t *nuclear_model,
        double *model_params,
        // answer
        double *matrix
)
{
    libgrpp_shell_t *shell_A = libgrpp_new_shell(origin_A, *L_A, *num_primitives_A, coeffs_A, alpha_A);
    libgrpp_shell_t *shell_B = libgrpp_new_shell(origin_B, *L_B, *num_primitives_B, coeffs_B, alpha_B);

    libgrpp_nuclear_attraction_integrals(shell_A, shell_B, charge_origin, *charge, *nuclear_model, model_params, matrix);

    libgrpp_delete_shell(shell_A);
    libgrpp_delete_shell(shell_B);
}


/*
 * Fortran interface to the nuclear models
 */

void libgrpp_estimate_nuclear_rms_radius_johnson_1985_(int32_t *A, double *R_rms)
{
    *R_rms = libgrpp_estimate_nuclear_rms_radius_johnson_1985(*A);
}


void libgrpp_estimate_nuclear_rms_radius_golovko_2008_(int32_t *A, double *R_rms)
{
    *R_rms = libgrpp_estimate_nuclear_rms_radius_golovko_2008(*A);
}


void libgrpp_estimate_fermi_model_parameters_(double *R_rms, double *c, double *a, int32_t *err_code)
{
    *err_code = (int32_t) libgrpp_estimate_fermi_model_parameters(*R_rms, c, a);
}


void libgrpp_charge_density_ball_(double *r, double *Z, double *R_rms, double *rho)
{
    *rho = libgrpp_charge_density_ball(*r, *Z, *R_rms);
}


void libgrpp_charge_density_gaussian_(double *r, double *Z, double *R_rms, double *rho)
{
    *rho = libgrpp_charge_density_gaussian(*r, *Z, *R_rms);
}


void libgrpp_charge_density_fermi_(double *r, double *Z, double *c, double *a, double *rho)
{
    *rho = libgrpp_charge_density_fermi(*r, *Z, *c, *a);
}


void libgrpp_charge_density_fermi_bubble_(double *r, double *Z, double *c, double *a, double *k, double *rho)
{
    *rho = libgrpp_charge_density_fermi_bubble(*r, *Z, *c, *a, *k);
}


void libgrpp_coulomb_potential_point_(double *r, double *Z, double *potential)
{
    *potential = libgrpp_coulomb_potential_point(*r, *Z);
}


void libgrpp_coulomb_potential_ball_(double *r, double *Z, double *R_rms, double *potential)
{
    *potential = libgrpp_coulomb_potential_ball(*r, *Z, *R_rms);
}


void libgrpp_coulomb_potential_gaussian_(double *r, double *Z, double *R_rms, double *potential)
{
    *potential = libgrpp_coulomb_potential_gaussian(*r, *Z, *R_rms);
}


void libgrpp_coulomb_potential_fermi_(double *r, double *Z, double *c, double *a, double *potential)
{
    *potential = libgrpp_coulomb_potential_fermi(*r, *Z, *c, *a);
}


void libgrpp_coulomb_potential_fermi_bubble_(double *r, double *Z, double *c, double *a, double *k, double *potential)
{
    *potential = libgrpp_coulomb_potential_fermi_bubble(*r, *Z, *c, *a, *k);
}


void libgrpp_rms_radius_fermi_(int32_t *Z, double *c, double *a, double *r_rms)
{
    *r_rms = libgrpp_rms_radius_fermi(*Z, *c, *a);
}


void libgrpp_rms_radius_fermi_bubble_(int32_t *Z, double *c, double *a, double *k, double *r_rms)
{
    *r_rms = libgrpp_rms_radius_fermi_bubble(*Z, *c, *a, *k);
}


