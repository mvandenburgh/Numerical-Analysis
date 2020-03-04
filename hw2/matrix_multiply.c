#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 1024

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
double randomNumber();


int main() {
    time_t startTime = time(NULL);
    srand( (unsigned) time(NULL) ); // use current time as seed for rng
    fillMatrix(A);
    fillMatrix(B);
    naiveMatrixMultiply(A,B,C);
    time_t runningTime = time(NULL) - startTime;
    printMatrix(A);
    printf("\n\n");
    printMatrix(B);
    printf("\n\n");
    printMatrix(C);
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

double randomNumber() {
    double random_value;

    random_value = ((double)rand()/(double)(RAND_MAX/2)) - 1.0;

    return random_value;
}

void naiveMatrixMultiply(double A[N][N], double B[N][N], double C[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double sum = 0;
            for (int k = 0; k < N; k++) {
                sum = sum + A[i][k] * B[k][j];
            }
            C[i][j] = sum;
        }
    }
}