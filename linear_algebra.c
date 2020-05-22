#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

/**
 * Add two nxn matrices together.
 * @param n dimemsion of nxn matrices
 */
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

/**
 * Subtract two nxn matrices
 * @param n dimension of nxn matrices
 */
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

/**
 * Multiplies two matrices of dimension N x N using naive O(N^3) algorithm
 */
void naive_matrix_multiply(int N, double** A, double** B, double** C) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C[i][j] = 0.0;
            for (int k = 0; k < N; k++) {
                C[i][j] = C[i][j] + A[i][k] * B[k][j];
            }
        }
    }
}

/**
 * Multiples two matrices of dimension N x N using Strassen's
 * algorithm (runs in O(N^2.8))
 * @todo fix gcc type warnings by switching to double** for matrices
 */
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

/**
 * Fills a matrix with random numbers.
 * @param matrix 2d array representing a matrix
 * @param rows number of rows in matrix
 * @param cols number of columns in matrix
 */
void fillMatrix(double** matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        matrix[i] = calloc(cols, sizeof(double));
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = (rand() % (9 - 0 + 1)) + 0;
        }
    }
}

/**
 * Free a matrix allocated by malloc(3).
 * @param matrix 2d array representing matrix
 * @param n dimension of nxn matrix to free
 */
void free_matrix(double** matrix, int n) {
    for (int i = 0; i < n; i++)
        free(matrix[i]);
    free(matrix);
}

/**
 * Perform LU Decomposition on a matrix of size n x n.
 * Caller must free both returned matrices after use.
 * @param matrix 2d array representing a matrix
 * @param n dimensions of nxn matrix
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
 * @param matrix 2d array representing a matrix
 * @param n dimensions of nxn matrix
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

/**
 * Calculates the determinant of an nxn matrix.
 * @param matrix 2d array representing a matrix
 * @param n dimensions of nxn matrix
 */
int det(double** matrix, int n) {
	if (n == 1) 
        return **matrix;
	else {
		int sign = 1;
		double ret = 0;
		double **sub;
		for (int i = 0; i < n; i++) {
			sub = malloc(sizeof(double *)*(n-1));
			for (int j = 0; j < n-1; j++) {
                sub[j] = malloc(sizeof(double)*(n-1));
            }
			int m = 0;
			for (int j = 0; j < n; j++) {
				if (j != i) {
					for (int k = 0; k < n-1; k++) {
                        sub[m][k] = matrix[j][k+1];
                    }
					m++;
				}
			}
			ret += matrix[i][0]*det(sub, n-1)*sign;
			sign = -sign;
			free(sub);
		}
		return ret;
	}
}
static const float convergence_threshold = 0.0000001f;

double** init_matrix(int n) {
    double** matrix = calloc(n, sizeof(double*));
    for (int i = 0; i < n; i++)
        matrix[i] = calloc(n, sizeof(double));
    return matrix;
}

void copy_matrix(double** dest, double** src, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dest[i][j] = src[i][j];
        }
    }
}

/**
 * Perform QR decomposition on a nxn matrix.
 * @param u matrix to perform qr decomposition on
 * @param n dimensions of nxn matrix
 * @param q allocated matrix to store Q part of decomposition
 * @param r allocated matrix to store R part of decomposition
 */
void qr_decomposition(double** u, int n, double** q, double** r) {
    long unsigned int i, j, k, ii;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            q[j][i] = u[j][i];
        }

        for (j = 0; j < i; j++) {
            double dot = 0;
            for (k = 0; k < n; k++) {
                dot += q[k][j] * u[k][i];
            }
            r[j][i] = dot;

            for (k = 0; k < n; k++) {
                q[k][i] = q[k][i] - r[j][i] * q[k][j];
            }
        }

        double l2_norm = 0;
        for (j = 0; j < n; j++)
            l2_norm += q[j][i] * q[j][i];
        r[i][i] = sqrtf(l2_norm);

        for (j = 0; j < n; j++) 
            q[j][i] = (1 / r[i][i]) * q[j][i];
    }
}

/**
 * Performs QR iterations on an N x N matrix, resulting in a new matrix
 * with the main diagonal containing the eigenvalues of the original.
 * @param matrix 2d array representing a matrix
 * @param result pre-allocated matrix to store result in
 * @param n dimension of nxn matrix
 * @param convergence_threshold the threshold at which a result is considered accurate enough (should be very small, for ex 0.0000001f)
 */
void qr_iterations(double** matrix, double** result, int n, float convergence_threshold) {
    long unsigned int i, j, k, ii;

    double** temp = init_matrix(n);
    double** q = init_matrix(n);
    double** r = init_matrix(n);

    copy_matrix(result, matrix, n);

    i = 0;
    while(1) {
        qr_decomposition(result, n, q, r);

        naive_matrix_multiply(n, r, q, temp);

        double total = 0.0;

        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                total += pow(temp[j][k] - result[j][k], 2);
            }
        }

        if (total < convergence_threshold || i > 1000000L)
            break;
        
        copy_matrix(result, temp, n);
        i++;
    }
    free_matrix(temp, n);
    free_matrix(q, n);
    free_matrix(r, n);
}

/**
 * Returns the approximate eigenvalues of a matrix using the
 * QR algorithm for numerically approximating eigenvalues.
 * @param matrix 2d array representing matrix
 * @param eigenvalues double array of length n to store the diagonal values in
 * @param n dimensions of nxn matrix / length of eigenvalues array
 */
void get_eigenvalues(double** matrix, double* eigenvalues, int n) {
    double** result = init_matrix(n);
    qr_iterations(matrix, result, n, 0.0000001f);
    for (int i = 0; i < n; i++) {
        eigenvalues[i] = result[i][i];
    }
    free_matrix(result, n);
}

void print_vector(double* vec, int n) {
    printf("[");
    for (int i = 0; i < n; i++) {
        printf("%f, ", vec[i]);
    }
    printf("]\n");
}


#define N 4 // size of matrix

int main() {
    time_t start_time = time(NULL);
    int rows = N;
    int cols = N;
    double** matrix = calloc(rows, sizeof(double*));
    fillMatrix(matrix, rows, cols);
    printf("Generated matrix:\n");
    print_matrix(matrix, rows, cols); printf("\n");
    double* eigenvalues = calloc(rows, sizeof(double));
    get_eigenvalues(matrix, eigenvalues, rows);

    printf("Eigenvalues:\n");
    print_vector(eigenvalues, rows);

    free(eigenvalues);
    double d = det(matrix, N);

    printf("\nDeterminant of matrix = %f\n", d);

    double*** lu_decomposition = lu_decomp(matrix, rows);
    double** lower = lu_decomposition[0];
    double** upper = lu_decomposition[1];
    double** inverted = inverse(matrix, rows);

    printf("\nInverse of matrix:\n");
    print_matrix(inverted, rows, cols);
    printf("\n\nRunning time: %lu second(s).\n", time(NULL) - start_time);

    // free allocated resources
    free_matrix(matrix, rows);
    free_matrix(inverted, rows);
    free_matrix(lower, rows);
    free_matrix(upper, rows);
    // free(lu_decomposition);
}