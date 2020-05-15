#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

double** matrixAddition(int n, double** A, double** B) {
    double** result = calloc(n, sizeof(double*));
    for (int i = 0; i < n; i++) {
        result[i] = calloc(n, sizeof(double));
        for (int j = 0; j < n; j++) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }
    return result;
}


double** matrixSubtraction(int n, double** A, double** B) {
    double** result = calloc(n, sizeof(double*));
    for (int i = 0; i < n; i++) {
        result[i] = calloc(n, sizeof(double));
        for (int j = 0; j < n; j++) {
            result[i][j] = A[i][j] - B[i][j];
        }
    }
    return result;
}


void naiveMatrixMultiply(int N, double** A, double** B, double** C) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C[i][j] = 0.0;
            for (int k = 0; k < N; k++) {
                C[i][j] = C[i][j] + A[i][k] * B[k][j];
            }
        }
    }
}

double** strassen(int n, double** A, double** B) {
    if (n <= 2) {
        double** result = calloc(2, sizeof(double*));
        result[0] = calloc(2, sizeof(double));
        result[1] = calloc(2, sizeof(double));
        result[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0];
        result[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1];
        result[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0];
        result[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1];
        
        return result;
    } 
    else {
        double (*a11)[n/2] = malloc(sizeof(double[n/2][n/2]));
        double (*a12)[n/2] = malloc(sizeof(double[n/2][n/2]));
        double (*a21)[n/2] = malloc(sizeof(double[n/2][n/2]));
        double (*a22)[n/2] = malloc(sizeof(double[n/2][n/2]));
        
        double (*b11)[n/2] = malloc(sizeof(double[n/2][n/2]));
        double (*b12)[n/2] = malloc(sizeof(double[n/2][n/2]));
        double (*b21)[n/2] = malloc(sizeof(double[n/2][n/2]));
        double (*b22)[n/2] = malloc(sizeof(double[n/2][n/2]));

        for (int i = 0; i < n/2; i++) {
            for (int j = 0; j < n/2; j++) {
                a11[i][j] = A[i][j];
                a12[i][j] = A[i][j + n/2];
                a21[i][j] = A[i + n/2][j];
                a22[i][j] = A[i + n/2][j + n/2];

                b11[i][j] = B[i][j];
                b12[i][j] = B[i][j + n/2];
                b21[i][j] = B[i + n/2][j];
                b22[i][j] = B[i + n/2][j + n/2];
            }
        }

        double (*temp1)[n];
        double (*temp2)[n];

        temp1 = matrixAddition(n/2, a11, a22);
        temp2 = matrixAddition(n/2, b11, b22);
        double (*m1)[n/2] = strassen(n/2, temp1, temp2);
        free(temp1); free(temp2);
        temp1 = matrixAddition(n/2, a21, a22);
        double (*m2)[n/2] = strassen(n/2, temp1, b11);
        free(temp1);
        temp1 = matrixSubtraction(n/2, b12, b22);
        double (*m3)[n/2] = strassen(n/2, a11, temp1);
        free(temp1);
        temp1 = matrixSubtraction(n/2, b21, b11);
        double (*m4)[n/2] = strassen(n/2, a22, temp1);
        free(temp1);
        temp1 = matrixAddition(n/2, a11, a12);
        double (*m5)[n/2] = strassen(n/2, temp1, b22);
        free(temp1);
        temp1 = matrixSubtraction(n/2, a21, a11);
        temp2 = matrixAddition(n/2, b11, b12);
        double (*m6)[n/2] = strassen(n/2, temp1, temp2);
        free(temp1); free(temp2);
        temp1 = matrixSubtraction(n/2, a12, a22);
        temp2 = matrixAddition(n/2, b21, b22);
        double (*m7)[n/2] = strassen(n/2, temp1, temp2);
        free(a11); free(a22); free(a21); free(a12); free(b11); free(b22); free(b21); free(b12); free(temp1); free(temp2);
        

        double (*c)[n] = malloc(sizeof(double[n][n]));

        for (int i = 0; i < n/2; i++) {
            for (int j = 0; j < n/2; j++) {
                c[i][j] = m1[i][j] + m4[i][j] - m5[i][j] + m7[i][j];
                c[i][j+n/2] = m3[i][j] + m5[i][j];
                c[i+n/2][j] = m2[i][j] + m4[i][j];
                c[i+n/2][j+n/2] = m1[i][j] - m2[i][j] + m3[i][j] + m6[i][j];
            }
        }
        free(m1); free(m2); free(m3); free(m4); free(m5); free(m6); free(m7);
        
        return (double **)c;
    }
}

void print_matrix(double** matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%.4f ", matrix[i][j]);
        }
        printf("\n");
    }
}

void fillMatrix(double** matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        matrix[i] = calloc(cols, sizeof(double));
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = (rand() % (9 - 0 + 1)) + 0;
        }
    }
}

void free_matrix(double** matrix, int rows) {
    for (int i = 0; i < rows; i++)
        free(matrix[i]);
    free(matrix);
}

/**
 * Perform LU Decomposition on a matrix of size n x n.
 * Caller must free both returned matrices after use.
 * @return 2 element array with 1st element being the lower triangular
 * matrix and the 2nd element being the upper triangular matrix.
 */
double*** lu_decomp(double** matrix, int n) {
    double** upper = calloc(n, sizeof(double));
    double** lower = calloc(n, sizeof(double));
    double sum;
    for (int i = 0; i < n; i++) {
        upper[i] = calloc(n, sizeof(double));
        lower[i] = calloc(n, sizeof(double));
    }
    for (int i = 0; i < n; i++) {
        for (int k = i; k < n; k++) {
            sum = 0.0;
            for (int j = 0; j < i; j++) {
                sum += (lower[i][j] * upper[j][k]);
            }
            upper[i][k] = matrix[i][k] - sum;
        }
        for (int k = i; k < n; k++) {
            if (i == k)
                lower[i][i] = 1;
            else {
                sum = 0.0;
                for (int j = 0; j < i; j++) {
                    sum += (lower[k][j] * upper[j][i]);
                }
                lower[k][i] = (matrix[k][i] - sum) / upper[i][i];
            }
        }
    }
    double*** ret = malloc(sizeof(lower) + sizeof(upper));
    ret[0] = lower;
    ret[1] = upper;
    return ret;
}

/**
 * Returns the inverse of the given matrix.
 * Caller must free the returned pointer after use.
 */
double** inverse(double** matrix, int n) {
    double** augmented = calloc(n, sizeof(double*));
    double temp;

    // Generate augmented matrix
    for (int i = 0; i < n; i++) {
        augmented[i] = calloc(2 * n, sizeof(double));
        for (int j = 0; j < n; j++) {
            augmented[i][j] = matrix[i][j];
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = n; j < 2 * n; j++) {
            if ((i + n) == j)
                augmented[i][j] = 1;
        }
    }

    // Perform Gauss-Jordan elimination
    for (int i = n - 1; i > 0; i--) {
        if (augmented[i-1][0] < augmented[i][0]) {
            double* tempp = augmented[i];
            augmented[i] = augmented[i-1];
            augmented[i-1] = tempp;
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (j != i) {
                temp = augmented[j][i] / augmented[i][i];
                for (int k = 0; k < 2 * n; k++) {
                    augmented[j][k] -= augmented[i][k] * temp;
                }
            }
        }
    }

    for (int i = 0; i < n; i++) {
        temp = augmented[i][i];
        for (int j = 0; j < 2 * n; j++) {
            augmented[i][j] = augmented[i][j] / temp;
        }
    }

    double** inverse = calloc(n, sizeof(double*));

    for (int i = 0; i < n; i++) {
        inverse[i] = calloc(n, sizeof(double));
        for (int j = 0; j < n; j++) {
            inverse[i][j] = augmented[i][j + n];
        }
    }
    free_matrix(augmented, n);
    return inverse;
}

int main() {
    int rows = 3;
    int cols = 3;
    double** matrix = calloc(rows, sizeof(double*));
    fillMatrix(matrix, rows, cols);
    double*** lu_decomposition = lu_decomp(matrix, rows);
    double** lower = lu_decomposition[0];
    double** upper = lu_decomposition[1];
    double** inverted = inverse(matrix, rows);


    print_matrix(matrix, rows, cols);
    printf("\n");
    print_matrix(inverted, rows, cols);

    // free allocated resources
    free_matrix(matrix, rows);
    free_matrix(inverted, rows);
    free_matrix(lower, rows);
    free_matrix(upper, rows);
    free(lu_decomposition);
}