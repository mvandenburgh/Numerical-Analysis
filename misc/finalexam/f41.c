#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

// Seeds for pseudo-RNG 
static long double seed0;
static long double seed1;

/**
 * Initialize seeds for pseudo-RNG
 */
void init_seeds(void) {
    seed0 = 0.15924135L;
    seed1 = 0.41214124L;
    // seed1 = (double)time(NULL);
    // while (seed1 >= 1) seed1 /= 10.0; // shift right a decimal place until > 1
}

/**
 * Generate a random floating point number between a and b
 * Algorithm sourced from Prof. Deng's lecture 2 notes
 */
long double random_float() {
    // long double temp = seed0;
    // seed0 = seed1;
    // if (seed1 + temp > 1.0)
    //     seed1 = fmodl((seed1 + temp) - 1.0, 1.0);
    // else
    //     seed1 = fmodl((seed1 + temp), 1.0);
    // return seed1;
    return (double)rand()/(double)(RAND_MAX);
}

double box_muller_transform(double mu, double sigma) {
    long double x1 = random_float();
    long double x2 = random_float();

    // printf("%lf ", x1);

    long double g1 = sqrt(-2 * log(x1)) * cos(2 * M_PI * x2);
    long double g2 = sqrt(-2 * log(x1)) * sin(2 * M_PI * x2);

    long double z1 = mu + g1 * sigma;
    long double z2 = mu + g2 * sigma;
    return z1;
}

/**
 * Initialize an empty matrix of size nxn
 */
double** init_matrix(int n) {
    double** matrix = calloc(n, sizeof(double*));
    for (int i = 0; i < n; i++)
        matrix[i] = calloc(n, sizeof(double));
    return matrix;
}

/**
 * Copy the contents of matrix src to another 
 * matrix dest.
 */
void copy_matrix(double** dest, double** src, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dest[i][j] = src[i][j];
        }
    }
}

/**
 * Frees an allocated matrix from memory.
 */
void free_matrix(double** matrix, int rows) {
    for (int i = 0; i < rows; i++)
        free(matrix[i]);
    free(matrix);
}

/**
 * Fills a matrix with appropriate values.
 */
void fillMatrix(double** matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        matrix[i] = calloc(cols, sizeof(double));
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = box_muller_transform(0, 1);
        }
    }
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

/**
 * Perform QR decomposition on a matrix, returning
 * the resulting matrices Q and R.
 */
void qr_decomposition(double** u, int n, double** q, double** r) {
    long unsigned int i, j, k, ii;

    for (i = 0; i < n; i++) {
        // q_i = u_i
        for (j = 0; j < n; j++) {
            q[j][i] = u[j][i];
        }

        for (j = 0; j < i; j++) {
        // r[j,i] = q_j^T * u_i 
            double dot = 0;
            for (k = 0; k < n; k++) {
                dot += q[k][j] * u[k][i];
            }
            r[j][i] = dot;

            // q_i = q_i - r[j,i] q_j
            for (k = 0; k < n; k++) {
                q[k][i] = q[k][i] - r[j][i] * q[k][j];
            }
        }

        // r[i,i] = |q_i|
        double l2_norm = 0;
        for (j = 0; j < n; j++)
            l2_norm += q[j][i] * q[j][i];
        r[i][i] = sqrtf(l2_norm);

        // q_i = q_i / r[i,i]
        for (j = 0; j < n; j++) 
            q[j][i] = (1 / r[i][i]) * q[j][i];
    }
}

/**
 * Perform the QR algorithm on a matrix, returning a resulting matrix
 * that has the eigenvalues of the original matrix as its diagnoal.
 */
void qr(double** matrix, double** result, int n) {
    long unsigned int i, j, k, ii;

    double** temp = init_matrix(n);
    double** q = init_matrix(n);
    double** r = init_matrix(n);

    copy_matrix(result, matrix, n);

    i = 0;
    while(1) {
        qr_decomposition(result, n, q, r);

        naiveMatrixMultiply(n, r, q, temp);

        double total = 0.0;

        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                total += pow(temp[j][k] - result[j][k], 2);
            }
        }

        if (total < 0.0000001f || i > 1000000L)
            break;
        
        copy_matrix(result, temp, n);
        i++;
    }
    free_matrix(temp, n);
    free_matrix(q, n);
    free_matrix(r, n);
}

/**
 * Returns the eigenvalues given a matrix that has 
 * been returned by the QR algorithm
 */
void get_eigenvalues(double** matrix, double* eigenvalues, int n) {
    double** result = init_matrix(n);
    qr(matrix, result, n);
    for (int i = 0; i < n; i++) {
        eigenvalues[i] = result[i][i];
    }
    free_matrix(result, n);
}

void print_eigenvalues(double* eigenvalues, int n) {
    printf("\nEigenvaules:\n");
    for (int i = 0; i < n; i++) {
        printf("%f ", eigenvalues[i]);
    }
    printf("\n");
}

void get_dominant_eigenvalue(double* eigenvalues, int n, double* dominant_eigenvalue) {
    int dominant = 0;
    for (int i = 1; i < n; i++) {
        if (abs(eigenvalues[i]) > abs(eigenvalues[dominant]))
            dominant = i;
    }
    FILE * temp = fopen("eigenvalues.txt", "w");
    for (int i=0; i < n; i++) {
        // for (int j = 0; j < n; j++) {
            fprintf(temp, "%lf ", eigenvalues[i]); //Write the data to a temporary file
        // }
        
    }
    fclose(temp);
    printf("Dominant eigenvalue = %f\n", eigenvalues[dominant]);
}

void print_matrix(double** matrix, int rows, int cols) {
    FILE * temp = fopen("matrix.txt", "w");
    fprintf(temp, "%d\n", rows);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%.4f ", matrix[i][j]);
            fprintf(temp, "%lf ", matrix[i][j]);
        }
        printf("\n");
        fprintf(temp, "\n");
    }
    fclose(temp);
}

int main() {
    int n = 12;
    double** matrix = init_matrix(n);
    fillMatrix(matrix, n, n);
    print_matrix(matrix, n,n);
    double* eigenvalues = calloc(n, sizeof(double));
    double* dominant;
    get_eigenvalues(matrix, eigenvalues, n);
    print_eigenvalues(eigenvalues, n);
    get_dominant_eigenvalue(eigenvalues, n, dominant);
    
}