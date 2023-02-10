/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2023 Alexander Oleynichenko
 */

#include "rpp.h"

#include <stdlib.h>
#include <string.h>

int skip_line(FILE *fp);
void print_rpp_expansion(FILE *out, libgrpp_potential_t *ecp);


libgrpp_grpp_t *read_grpp(char *path, int nuc_charge)
{
    FILE *inp_file;
    char buf[256];
    int z, n_blocks, n_arep;
    const int MAX_N_EXPONENTS = 100;
    const int MAX_N_CONTRACTED = 100;
    int pow_buf[MAX_N_EXPONENTS];
    double exp_buf[MAX_N_EXPONENTS];
    double coef_buf[MAX_N_CONTRACTED][MAX_N_EXPONENTS];

    // open file with basis set
    inp_file = fopen(path, "r");
    if (inp_file == NULL) {
        printf("Cannot open ecp file '%s\n'", path);
        return NULL;
    }

    // find the entry for the required element (nuclear charge)
    int found = 0;
    while (fscanf(inp_file, "%s", buf) == 1) {
        if (strcmp(buf, "*") == 0) {

            if (fscanf(inp_file, "%d%d%d", &z, &n_blocks, &n_arep) != 3) {
                printf("Error while reading ecp file for Z = %d\n", nuc_charge);
                return NULL;
            }

            if (z == nuc_charge) {
                found = 1;
                break;
            }
        }
    }

    if (found == 0) {
        return NULL;
    }

    // create template for the effective core potential
    libgrpp_grpp_t *grpp = libgrpp_new_grpp();

    int n_oc_shells = 0;
    libgrpp_shell_t *buf_oc_shells[MAX_NUM_OC_SHELLS];
    libgrpp_potential_t *buf_oc_potentials[MAX_NUM_OC_SHELLS];

    // read outercore pseudospinor expansions
    for (int iblock = 0; iblock < n_blocks; iblock++) {
        int n_alpha;
        int n_contr;
        double origin[] = {0, 0, 0};

        fscanf(inp_file, "%d%d", &n_alpha, &n_contr);
        skip_line(inp_file);

        for (int row = 0; row < n_alpha; row++) {
            fscanf(inp_file, "%lf", &exp_buf[row]);
            for (int col = 0; col < n_contr; col++) {
                fscanf(inp_file, "%lf", &coef_buf[col][row]);
            }
        }

        // add basis functions to the set
        for (int ifun = 0; ifun < n_contr; ifun++) {
            libgrpp_shell_t *shell = libgrpp_new_shell(origin, iblock, n_alpha, coef_buf[ifun], exp_buf);
            buf_oc_shells[n_oc_shells++] = shell;
        }
    }

    // read ECP expansions
    //grpp->n_oc_shells = 0;
    n_oc_shells = 0;
    for (int iblock = 0; iblock < n_arep; iblock++) {
        int n_prim;
        int n_oc;
        int L = iblock;

        fscanf(inp_file, "%d%d", &n_prim, &n_oc);
        skip_line(inp_file);


        for (int row = 0; row < n_prim; row++) {

            // power, exponent
            fscanf(inp_file, "%d%lf", &pow_buf[row], &exp_buf[row]);

            // arep, esop
            fscanf(inp_file, "%lf%lf", &coef_buf[0][row], &coef_buf[1][row]);

            // outercore potentials
            for (int ioc = 0; ioc < n_oc; ioc++) {
                fscanf(inp_file, "%lf", &coef_buf[2 + ioc][row]);
            }

            skip_line(inp_file);
        }

        // add potentials to the GRECP structure

        libgrpp_potential_t *arep = libgrpp_new_potential(L, 0, n_prim, pow_buf, coef_buf[0], exp_buf);
        if (L == n_arep - 1) { // special case of the U_L radially-local potential
            libgrpp_grpp_set_local_potential(grpp, arep);
        }
        else {
            libgrpp_grpp_add_averaged_potential(grpp, arep);
        }

        if (L == 0) {
            // no SO term for angular momentum S => put NULL
            libgrpp_grpp_add_spin_orbit_potential(grpp, NULL);
        }
        else {
            libgrpp_potential_t *esop = libgrpp_new_potential(L, 0, n_prim, pow_buf, coef_buf[1], exp_buf);
            libgrpp_grpp_add_spin_orbit_potential(grpp, esop);
        }

        // save OC potentials to buffer
        for (int ioc = 0; ioc < n_oc; ioc++) {
            int J = abs(2 * L + ((ioc % 2 == 0) ? -1 : +1));
            libgrpp_potential_t *oc_pot = libgrpp_new_potential(L, J, n_prim, pow_buf, coef_buf[2 + ioc], exp_buf);
            buf_oc_potentials[n_oc_shells++] = oc_pot;
        }
    }

    // add pairs (OC potential, OC shell) to the GRPP
    for (int ioc = 0; ioc < n_oc_shells; ioc++) {
        libgrpp_grpp_add_outercore_potential(grpp, buf_oc_potentials[ioc], buf_oc_shells[ioc]);
    }

    fclose(inp_file);

    return grpp;
}


void print_grpp(FILE *out, libgrpp_grpp_t *grpp)
{
    char ang_mom_labels[] = "SPDFGHIKLMNOPQ";

    printf("\toutercore pseudospinors:\n\n");

    for (int ioc = 0; ioc < grpp->n_oc_shells; ioc++) {
        libgrpp_shell_t *shell = grpp->oc_shells[ioc];
        char L_symbol = ang_mom_labels[shell->L];

        for (int i = 0; i < shell->num_primitives; i++) {
            fprintf(out, "  %c%18.8e%20.8e\n", (i == 0) ? L_symbol : ' ', shell->alpha[i], shell->coeffs[i]);
        }

        printf("\n");
    }

    printf("\t\tradially-local potential:\n\n");
    print_rpp_expansion(stdout, grpp->U_L);
    printf("\n");

    printf("\t\tsemi-local averaged potentials:\n\n");
    for (int L = 0; L < grpp->n_arep; L++) {
        print_rpp_expansion(stdout, grpp->U_arep[L]);
        printf("\n");
    }

    printf("\t\tsemi-local spin-orbit potentials:\n\n");
    for (int L = 1; L < grpp->n_esop; L++) {
        print_rpp_expansion(stdout, grpp->U_esop[L]);
        printf("\n");
    }

    printf("\t\tnon-local outercore potentials:\n\n");
    for (int ioc = 0; ioc < grpp->n_oc_shells; ioc++) {
        print_rpp_expansion(stdout, grpp->U_oc[ioc]);
        printf("\n");
    }
}


void print_rpp_expansion(FILE *out, libgrpp_potential_t *ecp)
{
    char ang_mom_labels[] = "SPDFGHIKLMNOPQ";
    char L_symbol = ang_mom_labels[ecp->L];
    char J_symbol[8];

    if (ecp->J > 0) {
        sprintf(J_symbol, "%2d/2", ecp->J);
    }
    else {
        sprintf(J_symbol, "    ");
    }

    for (int i = 0; i < ecp->num_primitives; i++) {
        fprintf(out, "  %c%4s%4d%20.8e%20.8e\n",
                (i == 0) ? L_symbol : ' ',
                (i == 0) ? J_symbol : "",
                ecp->powers[i], ecp->alpha[i], ecp->coeffs[i]);
    }
}


int skip_line(FILE *fp)
{
    int c;

    while (c = fgetc(fp), c != '\n' && c != EOF);

    return c;
}
