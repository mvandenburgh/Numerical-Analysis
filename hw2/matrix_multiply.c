/**
 * C program that implements both the naive method and Strassen's method
 * for matrix multiplication. Tested on gcc version 7.4.0.
 * @author Michael VanDenburgh
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define N 4096 // dimensions of two randomly generated matrices
#define PRINT 1 // set to 0 to disable printing out the matrices
#define NAIVE 1 // set to 0 to disable naive method
#define STRASSEN 1 // set to 0 to disable strassen method

// N = 100, less than a second
// N = 500, less than a second
// N = 600, about 1 second
// N = 700, about 1 second
// N = 800, about 3 seconds
// N = 900, about 4 seconds
// N = 1024, about 14 seconds


// Strassen N = 4096, 372 seconds

static void printMatrixNxN(int n, double m[n][n]);
static void fillMatrix(double matrix[N][N]);
static int matrixCmp(int n, double m1[n][n], double m2[n][n]);
static void naiveMatrixMultiply(double A[N][N], double B[N][N], double C[N][N]);
static double** matrixAddition(int n, double A[n][n], double B[n][n]);
static double** matrixSubtraction(int n, double A[n][n], double B[n][n]);
static double** strassen(int n, double A[n][n], double B[n][n]);
static double randomNumber();

static double A[N][N];
static double B[N][N];
static double C[N][N];


int main(void) {
    time_t startTime;
    time_t runningTime;

    srand( (unsigned) time(NULL) ); // use current time as seed for rng


    fillMatrix(A);
    fillMatrix(B);
    #if PRINT
    printf("Generated matrix A:");
    printMatrixNxN(N, A);
    printf("\n\nGenerated matrix B:");
    printMatrixNxN(N, B);
    #endif

    #if NAIVE
    printf("\nRunning naive method...");
    startTime = time(NULL);
    naiveMatrixMultiply(A, B, C);
    runningTime = time(NULL) - startTime;
    // #if PRINT
    printf("\nResult:");
    printMatrixNxN(N, C);
    // #endif
    printf("Running time of naive method: %lu seconds\n\n", runningTime);
    #endif
    
    #if STRASSEN
    printf("\nRunning Strassen's method...");
    startTime = time(NULL);
    double (*s)[N] = strassen(N, A, B);
    runningTime = time(NULL) - startTime;
    // #if PRINT
    printf("\nResult:");
    printMatrixNxN(N, s);
    // #endif
    printf("Running time of Strassen method: %lu seconds\n\n", runningTime);
    free(s);
    #endif

    // if (matrixCmp(N, s, C)) printf("\nPASS: MATRICES ARE EQUAL.\n\n");
    // else printf("FAIL: MATRICES ARE NOT EQUAL.\n\n");
}



void fillMatrix(double matrix[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i][j] = randomNumber();
        }
    }
}

// TODO: code my own pseudo rng
double randomNumber() {
    double random_value;

    random_value = ((double)rand()/(double)(RAND_MAX/2)) - 1.0;

    return random_value;
}


//TODO: take in preallocated block of memory instead of malloc'ing inside function
double** matrixAddition(int n, double A[n][n], double B[n][n]) {
    double (*result)[n] = malloc(sizeof(double[n][n]));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }
    return result;
}

//TODO: take in preallocated block of memory instead of malloc'ing inside function
double** matrixSubtraction(int n, double A[n][n], double B[n][n]) {
    double (*result)[n] = malloc(sizeof(double[n][n]));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i][j] = A[i][j] - B[i][j];
        }
    }
    return result;
}


void naiveMatrixMultiply(double A[N][N], double B[N][N], double C[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C[i][j] = 0.0;
            for (int k = 0; k < N; k++) {
                C[i][j] = C[i][j] + A[i][k] * B[k][j];
            }
        }
    }
}

double** strassen(int n, double A[n][n], double B[n][n]) {
    if (n <= 2) {
        double (*result)[2] = malloc(sizeof(double[2][2]));
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


void printMatrixNxN(int n, double m[n][n]) {
    printf("\n---------------------------------------------------------------------------------\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", m[i][j]);
            fflush(stdout);
        }
        printf("\n");
    }
    printf("---------------------------------------------------------------------------------\n");
}

int matrixCmp(int n, double m1[n][n], double m2[n][n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (fabs(m1[i][j] - m2[i][j]) > 0.001) {
                printf("\n%f != %f @ [%d][%d]: ", m1[i][j], m2[i][j], i, j);
                return 0;
            }
        }
    }
    return 1;
}