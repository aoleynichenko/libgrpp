/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2022 Alexander Oleynichenko
 */

#include "utils.h"

#include <math.h>
#include <stdlib.h>

inline int int_max2(int x, int y)
{
    return (x > y) ? x : y;
}


inline int int_max3(int x, int y, int z)
{
    return int_max2(int_max2(x, y), z);
}


double *alloc_zeros_1d(int n)
{
    return (double *) calloc(n, sizeof(double));
}

double **alloc_zeros_2d(int n, int m)
{
    double **array = (double **) calloc(n, sizeof(double *));

    for (int i = 0; i < n; i++) {
        array[i] = (double *) calloc(m, sizeof(double));
    }

    return array;
}

void free_2d(double **array, int n)
{
    for (int i = 0; i < n; i++) {
        free(array[i]);
    }
    free(array);
}

/*
 * constant times a vector plus a vector:
 * y = a * x + y
 */
void libgrpp_daxpy(int n, double a, double *x, double *y)
{
    for (int i = 0; i < n; i++) {
        y[i] += a * x[i];
    }
}


/*
 * naive matrix multiplication
 */
void libgrpp_multiply_matrices(int M, int N, int K, double *A, double *B, double *C)
{
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            double sum = 0.0;
            for (int k = 0; k < K; k++) {
                sum += A[i * K + k] * B[k * N + j];
            }
            C[i * N + j] += sum;
        }
    }
}


double distance_squared(double *A, double *B)
{
    double dx = A[0] - B[0];
    double dy = A[1] - B[1];
    double dz = A[2] - B[2];
    return dx * dx + dy * dy + dz * dz;
}


double distance(double *A, double *B)
{
    return sqrt(distance_squared(A, B));
}
