#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 5


static double A[N][N];
static double B[N][N];
static double C[N][N];

void printMatrix(double matrix[N][N]);
void fillMatrix(double matrix[N][N]);
void naiveMatrixMultiply(double A[N][N], double B[N][N], double C[N][N]);
double randomNumber();


int main() {
    srand( (unsigned) time(NULL) ); // use current time as seed for rng
    fillMatrix(A);
    fillMatrix(B);
    naiveMatrixMultiply(A,B,C);
    // printMatrix(A);
    // printf("\n");
    // printMatrix(B);
    // printf("\n");
    printMatrix(C);
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