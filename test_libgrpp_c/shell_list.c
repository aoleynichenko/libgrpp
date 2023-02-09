/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2023 Alexander Oleynichenko
 */

#include "shell_list.h"

#include <stdio.h>
#include <stdlib.h>


void get_cartesian_label(char *buf, int n, int l, int m);


libgrpp_shell_t **construct_shell_list(molecule_t *mol, basis_set_t **basis_sets, int *num_shells)
{
    // pre-calculate number of shells
    *num_shells = 0;
    for (int iatom = 0; iatom < mol->n_atoms; iatom++) {
        int nuc_charge = mol->charges[iatom];
        basis_set_t *bs = basis_sets[nuc_charge];
        *num_shells = *num_shells + bs->n_fun;
    }

    libgrpp_shell_t **shell_list = (libgrpp_shell_t **) malloc(sizeof(libgrpp_shell_t *) * (*num_shells));

    // generate atom-centered shells from basis functions
    int shell_count = 0;
    for (int iatom = 0; iatom < mol->n_atoms; iatom++) {
        int nuc_charge = mol->charges[iatom];

        double origin[3];
        origin[0] = mol->coord_x[iatom];
        origin[1] = mol->coord_y[iatom];
        origin[2] = mol->coord_z[iatom];

        basis_set_t *bs = basis_sets[nuc_charge];

        for (int ifun = 0; ifun < bs->n_fun; ifun++) {
            basis_fun_t *fun = bs->fun_list[ifun];
            shell_list[shell_count] = libgrpp_new_shell(origin, fun->L, fun->n_primitives, fun->coeffs, fun->alpha);
            shell_count++;
        }
    }

    return shell_list;
}


void delete_shell_list(libgrpp_shell_t **shell_list, int num_shells)
{
    for (int ishell = 0; ishell < num_shells; ishell++) {
        libgrpp_delete_shell(shell_list[ishell]);
    }
    free(shell_list);
}


int calculate_basis_dim(libgrpp_shell_t **shell_list, int num_shells)
{
    int basis_dim = 0;

    for (int ishell = 0; ishell < num_shells; ishell++) {
        libgrpp_shell_t *shell = shell_list[ishell];
        basis_dim += libgrpp_get_shell_size(shell);
    }

    return basis_dim;
}


void print_shell_list(libgrpp_shell_t **shell_list, int num_shells)
{
    int count = 1;
    char cart_label[32];

    printf("\n");
    printf("\t\tlist of basis shells:\n");
    printf("\t\t---------------------\n\n");
    printf("    \t           origin\n");
    printf("   #\t     x       y       z  \tcartesian\n");
    for (int ishell = 0; ishell < num_shells; ishell++) {
        libgrpp_shell_t *shell = shell_list[ishell];
        int shell_size = libgrpp_get_shell_size(shell);
        for (int icart = 0; icart < shell_size; icart++) {
            int n = shell->cart_list[3 * icart + 0];
            int l = shell->cart_list[3 * icart + 1];
            int m = shell->cart_list[3 * icart + 2];
            get_cartesian_label(cart_label, n, l, m);
            printf("%4d\t%8.3f%8.3f%8.3f\t%s\n", count++, shell->origin[0], shell->origin[1], shell->origin[2], cart_label);
        }
    }
    printf("\n");
}


void get_cartesian_label(char *buf, int n, int l, int m)
{
    int count = 0;

    for (int i = 0; i < n; i++) {
        buf[count++] = 'x';
    }
    for (int i = 0; i < l; i++) {
        buf[count++] = 'y';
    }
    for (int i = 0; i < m; i++) {
        buf[count++] = 'z';
    }

    // special case: S function
    if (count == 0) {
        buf[count++] = 's';
    }

    buf[count] = '\0';
}
