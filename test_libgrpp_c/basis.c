/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2022 Alexander Oleynichenko
 */


#include "basis.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int skip_line(FILE *fp);

void radial_gto_norm_factor_(int *L, int *num_primitives, double *coeffs, double *alpha, double *norm);


basis_fun_t *new_basis_function(int L, int n_primitives, double *coeffs, double *alpha)
{
    basis_fun_t *fun = (basis_fun_t *) malloc(sizeof(basis_fun_t) * 1);

    fun->L = L;
    fun->n_primitives = n_primitives;
    fun->coeffs = (double *) malloc(sizeof(double) * n_primitives);
    fun->alpha = (double *) malloc(sizeof(double) * n_primitives);
    for (int i = 0; i < n_primitives; i++) {
        fun->coeffs[i] = coeffs[i];
        fun->alpha[i] = alpha[i];
    }

    return fun;
}


void delete_basis_function(basis_fun_t *fun)
{
    free(fun->alpha);
    free(fun->coeffs);
    free(fun);
}


/*
typedef struct {
    int n_fun;
    basis_fun_t *fun_list;
} basis_set_t;
*/

basis_set_t *new_basis_set()
{
    basis_set_t *set = (basis_set_t *) malloc(sizeof(basis_set_t) * 1);

    set->n_fun = 0;
    set->capacity = 1;
    set->fun_list = (basis_fun_t **) malloc(sizeof(basis_fun_t *) * set->capacity);

    return set;
}


void delete_basis_set(basis_set_t *set)
{
    for (int i = 0; i < set->n_fun; i++) {
        delete_basis_function(set->fun_list[i]);
    }
    free(set->fun_list);
    free(set);
}


void add_basis_function(basis_set_t *set, basis_fun_t *fun)
{
    // check if the memory must be allocated
    if (set->n_fun == set->capacity) {
        set->fun_list = (basis_fun_t **) realloc(set->fun_list, 2 * set->capacity * sizeof(basis_fun_t *));
        set->capacity *= 2;
    }

    set->fun_list[set->n_fun] = fun;
    set->n_fun++;
}


basis_set_t *read_basis_set(char *path, int nuc_charge)
{
    FILE *inp_file;
    char buf[256];
    int z, n_blocks;
    const int MAX_N_EXPONENTS = 100;
    const int MAX_N_CONTRACTED = 100;
    double exp_buf[MAX_N_EXPONENTS];
    double coef_buf[MAX_N_CONTRACTED][MAX_N_EXPONENTS];

    // open file with basis set
    inp_file = fopen(path, "r");
    if (inp_file == NULL) {
        printf("Cannot open basis file '%s\n'", path);
        return NULL;
    }

    // find the entry for the required element (nuclear charge)
    int found = 0;
    while (fscanf(inp_file, "%s", buf) == 1) {
        if (strcmp(buf, "*") == 0) {

            if (fscanf(inp_file, "%d%d", &z, &n_blocks) != 2) {
                printf("Error while reading basis file for Z = %d\n", nuc_charge);
                return NULL;
            }

            if (z == nuc_charge) {
                found = 1;
                break;
            }
        }
    }

    if (found == 0) {
        printf("Error while reading basis file\n");
        printf("Basis set for Z = %d not found\n", nuc_charge);
        return NULL;
    }

    // create template for the basis set
    basis_set_t *set = new_basis_set();

    // read the entry block-wise
    for (int iblock = 0; iblock < n_blocks; iblock++) {
        int n_primitives;
        int n_contracted;

        fscanf(inp_file, "%d%d", &n_primitives, &n_contracted);

        for (int row = 0; row < n_primitives; row++) {
            fscanf(inp_file, "%lf", &exp_buf[row]);
            for (int col = 0; col < n_contracted; col++) {
                fscanf(inp_file, "%lf", &coef_buf[col][row]);
            }
            skip_line(inp_file);
        }

        // add basis functions to the set
        // each function must be normalized to unity
        for (int ifun = 0; ifun < n_contracted; ifun++) {

            double norm_factor;
            radial_gto_norm_factor_(&iblock, &n_primitives, coef_buf[ifun], exp_buf, &norm_factor);
            for (int i = 0; i < n_primitives; i++) {
                coef_buf[ifun][i] *= norm_factor;
            }

            add_basis_function(set, new_basis_function(iblock, n_primitives, coef_buf[ifun], exp_buf));
        }
    }

    fclose(inp_file);

    return set;
}


void print_basis_set(FILE *out, basis_set_t *set)
{
    double const PRINT_THRESH = 1e-8;
    char ang_mom_labels[] = "SPDFGHIKLMNOPQ";

    for (int ifun = 0; ifun < set->n_fun; ifun++) {
        basis_fun_t *fun = set->fun_list[ifun];
        char L_symbol = ang_mom_labels[fun->L];
        int n_zero = 0;

        for (int i = 0; i < fun->n_primitives; i++) {
            if (fabs(fun->coeffs[i]) < PRINT_THRESH) {
                n_zero++;
                continue;
            }
            fprintf(out, "  %c%18.8e%20.8e\n", (i - n_zero == 0) ? L_symbol : ' ', fun->alpha[i], fun->coeffs[i]);
        }

        printf("\n");
    }
}

