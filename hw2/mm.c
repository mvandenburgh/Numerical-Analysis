#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 16

// N = 100, less than a second
// N = 500, less than a second
// N = 600, about 1 second
// N = 700, about 1 second
// N = 800, about 3 seconds
// N = 900, about 4 seconds
// N = 1024, about 14 seconds

static double A[N][N];
static double B[N][N];
static double C[N][N];

void printMatrix(double matrix[N][N]);
void fillMatrix(double matrix[N][N]);
void naiveMatrixMultiply(double A[N][N], double B[N][N], double C[N][N]);
double** matrixAddition(int n, double A[n][n], double B[n][n]);
double** matrixSubtraction(int n, double A[n][n], double B[n][n]);
// double** strassen(double** A, double** B, int n);

double randomNumber();

void printMatrixNxN(int n, double m[n][n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", m[i][j]);
            fflush(stdout);
        }
        printf("\n");
    }
}

int areEqualMatrices(int n, double m1[n][n], double m2[n][n]) {
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


double** strassen(int n, double A[n][n], double B[n][n]) {
    if (n <= 2) {
        double m1 = (A[0][0] + A[1][1]) * (B[0][0] + B[1][1]);
        double m2 = (A[1][0] + A[1][1]) *  B[0][0];
        double m3 = A[0][0] * (B[0][1] - B[1][1]);
        double m4 = A[1][1] * (B[1][0] - B[0][0]);
        double m5 = (A[0][0] + A[0][1]) *  B[1][1];
        double m6 = (A[1][0] - A[0][0]) * (B[0][0] + B[0][1]);
        double m7 = (A[0][1] - A[1][1]) * (B[1][0] + B[1][1]);
        
        double (*result)[2] = malloc(sizeof(double[2][2]));
        result[0][0] = m1 + m4 - m5 + m7;
        result[0][1] = m3 + m5;
        result[1][0] = m2 + m4;
        result[1][1] = m1 - m2 + m3 + m6;
        // result[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0];
        // result[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1];
        // result[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0];
        // result[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1];
        
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

        
        // printMatrix(A);

        double (*temp1)[n];
        double (*temp2)[n];

        temp1 = matrixAddition(n/2, a11, a22);
        temp2 = matrixAddition(n/2, b11, b22);
        double (*d1)[n/2] = strassen(n/2, temp1, temp2);
        free(temp1); free(temp2);
        temp1 = matrixAddition(n/2, a21, a22);
        double (*d2)[n/2] = strassen(n/2, temp1, b11);
        free(temp1);
        temp1 = matrixSubtraction(n/2, b12, b22);
        double (*d3)[n/2] = strassen(n/2, a11, temp1);
        free(temp1);
        temp1 = matrixSubtraction(n/2, b21, b11);
        double (*d4)[n/2] = strassen(n/2, a22, temp1);
        free(temp1);
        temp1 = matrixAddition(n/2, a11, a12);
        double (*d5)[n/2] = strassen(n/2, temp1, b22);
        free(temp1);
        temp1 = matrixSubtraction(n/2, a21, a11);
        temp2 = matrixAddition(n/2, b11, b12);
        double (*d6)[n/2] = strassen(n/2, temp1, temp2);
        free(temp1); free(temp2);
        temp1 = matrixSubtraction(n/2, a12, a22);
        temp2 = matrixAddition(n/2, b21, b22);
        double (*d7)[n/2] = strassen(n/2, temp1, temp2);
        free(a11); free(a22); free(a21); free(a12); free(b11); free(b22); free(b21); free(b12); free(temp1); free(temp2);
        
        
        
        double (*result)[n] = malloc(sizeof(double[n][n]));

        for (int i = 0; i < n/2; i++) {
            for (int j = 0; j < n/2; j++) {
                result[i][j] = d1[i][j] + d4[i][j] - d5[i][j] + d7[i][j];
                result[i][j+n/2] = d3[i][j] + d5[i][j];
                result[i+n/2][j] = d2[i][j] + d4[i][j];
                result[i+n/2][j+n/2] = d1[i][j] - d2[i][j] + d3[i][j] + d6[i][j];
            }
        }
        free(d1); free(d2); free(d3); free(d4); free(d5); free(d6); free(d7);
        return result;
    }
}

int main() {
    time_t startTime = time(NULL);

    srand( (unsigned) time(NULL) ); // use current time as seed for rng

    fillMatrix(A);
    fillMatrix(B);
    
    naiveMatrixMultiply(A, B, C);
    double (*s)[N] = strassen(N, A, B);

    time_t runningTime = time(NULL) - startTime;

    printf("Naive:\n\n");
    printMatrixNxN(N, C);
    
    printf("\n");
    
    printf("Strassen:\n\n");
    printMatrixNxN(N, s);
    free(s);
    
    // if (areEqualMatrices(N, s, C)) printf("\nPASS: MATRICES ARE EQUAL.\n\n");
    // else printf("FAIL: MATRICES ARE NOT EQUAL.\n\n");

    printf("\nRunning time of program: %lu seconds\n", runningTime);

}


void printMatrix(double matrix[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
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

