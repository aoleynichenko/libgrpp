#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int in_range(int num, int min, int max);

int read_matrix(char *path, int *dim, double **matrix);

void print_usage(char *path_to_binary);

int compare_matrices(int dim, double *matrix1, double *matrix2, double threshold, int *diff_i, int *diff_j);


int main(int argc, char **argv)
{
    /*
     * parse command line arguments
     */
    if (argc < 3 || argc > 4) {
        print_usage(argv[0]);
        return 1;
    }

    double threshold = 1e-10;
    char *path_matrix_1 = argv[1];
    char *path_matrix_2 = argv[2];

    if (argc == 4) {
        int n_read = sscanf(argv[1], "%lf", &threshold);
        if (n_read != 1) {
            printf("error while parsing the threshold value\n");
            return 1;
        }

        path_matrix_1 = argv[2];
        path_matrix_2 = argv[3];
    }

    /*
     * read matrices to be compared
     */
    int dim1;
    double *matrix1 = NULL;
    if (read_matrix(path_matrix_1, &dim1, &matrix1) == -1) {
        printf("error while reading the first matrix\n");
        return 1;
    }

    int dim2;
    double *matrix2 = NULL;
    if (read_matrix(path_matrix_2, &dim2, &matrix2) == -1) {
        printf("error while reading the second matrix\n");
        return 1;
    }

    if (dim1 != dim2) {
        printf("dimensions of matrices do not coincide (%d vs %d)\n", dim1, dim2);
        return 1;
    }

    /*
     * compare matrices
     */
    int diff_i, diff_j;

    if (compare_matrices(dim1, matrix1, matrix2, threshold, &diff_i, &diff_j) == 1) {
        printf("matrices are equal\n");
        return 0;
    }
    else {
        double z1 = matrix1[diff_i * dim1 + diff_j];
        double z2 = matrix2[diff_i * dim1 + diff_j];
        printf("matrices are different:\n");
        printf("matrix1 [%d][%d] %16.12f\n", diff_i + 1, diff_j + 1, z1);
        printf("matrix2 [%d][%d] %16.12f\n", diff_i + 1, diff_j + 1, z2);
        return 1;
    }
}


void print_usage(char *path_to_binary)
{
    printf("Usage: %s [<double thresh>] <matrix_file_1> <matrix_file_2>\n", path_to_binary);
}


int read_matrix(char *path, int *dim, double **matrix)
{
    FILE *f = fopen(path, "r");
    if (f == NULL) {
        return -1;
    }

    // read dimension of a matrix
    if (fscanf(f, "%d", dim) != 1) {
        return -1;
    }

    // allocate memory
    *matrix = (double *) calloc((*dim) * (*dim), sizeof(double));
    if (*matrix == NULL) {
        return -1;
    }

    // read matrix
    int i, j;
    double val;
    while (fscanf(f, "%d%d%lf", &i, &j, &val) == 3) {
        if (!in_range(i, 1, *dim) || !in_range(j, 1, *dim)) {
            return -1;
        }

        (*matrix)[(*dim) * (i - 1) + (j - 1)] = val;
    }

    fclose(f);

    return 0;
}


/**
 * compares two complex matrices element-by-element.
 *
 * if matrices are not equal to each other, the position
 * of the difference is returned via the 'diff_pos' argument.
 */
int compare_matrices(int dim, double *matrix1, double *matrix2, double threshold, int *diff_i, int *diff_j)
{
    *diff_i = 0;
    *diff_j = 0;

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            double z1 = matrix1[i * dim + j];
            double z2 = matrix2[i * dim + j];

            if (fabs(z1 - z2) > threshold) {
                *diff_i = i;
                *diff_j = j;
                return 0;
            }
        }
    }

    return 1;
}


/**
 * checks if a number is in a range [min,max]
 */
int in_range(int num, int min, int max)
{
    if (min <= num && num <= max) {
        return 1;
    }
    return 0;
}
