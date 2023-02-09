/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2023 Alexander Oleynichenko
 */

#include "elements.h"

#include <assert.h>
#include <ctype.h>
#include <string.h>

/*
 * mass numbers are given for the most abundant isotopes
 * (or for the most long-lived for radioactive elements)
 */

element_t periodic_table[] = {
        {0,   "Gh", 0},
        {1,   "H",  1},
        {2,   "He", 4},
        {3,   "Li", 7},
        {4,   "Be", 9},
        {5,   "B",  11},
        {6,   "C",  12},
        {7,   "N",  14},
        {8,   "O",  16},
        {9,   "F",  19},
        {10,  "Ne", 20},
        {11,  "Na", 23},
        {12,  "Mg", 24},
        {13,  "Al", 27},
        {14,  "Si", 28},
        {15,  "P",  31},
        {16,  "S",  32},
        {17,  "Cl", 35},
        {18,  "Ar", 40},
        {19,  "K",  39},
        {20,  "Ca", 40},
        {21,  "Sc", 45},
        {22,  "Ti", 48},
        {23,  "V",  51},
        {24,  "Cr", 52},
        {25,  "Mn", 55},
        {26,  "Fe", 56},
        {27,  "Co", 59},
        {28,  "Ni", 58},
        {29,  "Cu", 63},
        {30,  "Zn", 64},
        {31,  "Ga", 69},
        {32,  "Ge", 74},
        {33,  "As", 75},
        {34,  "Se", 80},
        {35,  "Br", 79},
        {36,  "Kr", 84},
        {37,  "Rb", 85},
        {38,  "Sr", 88},
        {39,  "Y",  89},
        {40,  "Zr", 90},
        {41,  "Nb", 93},
        {42,  "Mo", 98},
        {43,  "Tc", 98},
        {44,  "Ru", 102},
        {45,  "Rh", 100},
        {46,  "Pd", 106},
        {47,  "Ag", 107},
        {48,  "Cd", 114},
        {49,  "In", 115},
        {50,  "Sn", 120},
        {51,  "Sb", 121},
        {52,  "Te", 130},
        {53,  "I",  127},
        {54,  "Xe", 132},
        {55,  "Cs", 133},
        {56,  "Ba", 138},
        {57,  "La", 139},
        {58,  "Ce", 140},
        {59,  "Pr", 141},
        {60,  "Nd", 142},
        {61,  "Pm", 145},
        {62,  "Sm", 152},
        {63,  "Eu", 153},
        {64,  "Gd", 158},
        {65,  "Tb", 159},
        {66,  "Dy", 164},
        {67,  "Ho", 165},
        {68,  "Er", 166},
        {69,  "Tm", 169},
        {70,  "Yb", 174},
        {71,  "Lu", 175},
        {72,  "Hf", 180},
        {73,  "Ta", 181},
        {74,  "W",  184},
        {75,  "Re", 187},
        {76,  "Os", 192},
        {77,  "Ir", 193},
        {78,  "Pt", 195},
        {79,  "Au", 197},
        {80,  "Hg", 202},
        {81,  "Tl", 205},
        {82,  "Pb", 208},
        {83,  "Bi", 209},
        {84,  "Po", 209},
        {85,  "At", 210},
        {86,  "Rn", 222},
        {87,  "Fr", 223},
        {88,  "Ra", 226},
        {89,  "Ac", 227},
        {90,  "Th", 232},
        {91,  "Pa", 231},
        {92,  "U",  238},
        {93,  "Np", 237},
        {94,  "Pu", 244},
        {95,  "Am", 243},
        {96,  "Cm", 247},
        {97,  "Bk", 247},
        {98,  "Cf", 251},
        {99,  "Es", 252},
        {100, "Fm", 257},
        {101, "Md", 258},
        {102, "Nb", 259},
        {103, "Lr", 266},
        {104, "Rf", 267},
        {105, "Db", 268},
        {106, "Sg", 269},
        {107, "Bh", 270},
        {108, "Hs", 277},
        {109, "Mt", 278},
        {110, "Ds", 281},
        {111, "Rg", 282},
        {112, "Cn", 285},
        {113, "Nh", 286},
        {114, "Fl", 289},
        {115, "Mc", 290},
        {116, "Lv", 293},
        {117, "Ts", 294},
        {118, "Og", 294},
        {119, "E119", 119+184},
        {120, "E120", 120+184},
        {121, "E121", 121+184},
        {122, "E122", 122+184},
        {123, "E123", 123+184},
        {124, "E124", 124+184},
        {125, "E125", 125+184},
        {126, "E126", 126+184},
        {127, "E127", 127+184},
        {128, "E128", 128+184},
        {129, "E129", 129+184},
        {130, "E130", 130+184}
};


void get_element_symbol(int z, char *sym)
{
    assert(z >= 0 && z <= N_CHEM_ELEMENTS);

    strcpy(sym, periodic_table[z].sym);
}


int get_element_nuc_charge(char *sym)
{
    const int n_elements = sizeof(periodic_table) / sizeof(element_t);
    char buf[MAX_ELEMENT_SYMBOL];

    // cast 'sym' to the 'standard form' (only the first letter is capitalized)
    strncpy(buf, sym, MAX_ELEMENT_SYMBOL);
    buf[MAX_ELEMENT_SYMBOL - 1] = '\0';
    buf[0] = toupper(buf[0]);
    for (int i = 1; i < MAX_ELEMENT_SYMBOL - 1; i++) {
        buf[i] = tolower(buf[i]);
    }

    for (int i = 0; i < n_elements; i++) {
        if (strcmp(periodic_table[i].sym, buf) == 0) {
            return i;
        }
    }

    return -1;
}


int get_element_mass_number_abundant(int z)
{
    return periodic_table[z].mass_number;
}

