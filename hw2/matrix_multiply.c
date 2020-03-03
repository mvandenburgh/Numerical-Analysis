#include <stdio.h>
#include <stdlib.h>

#define N 100


static double matrix[N][N];

void printMatrix(double matrix[N][N]);
void fillMatrix(double matrix[N][N]);
double randomNumber();

int main() {
    fillMatrix(matrix);
    printMatrix(matrix);
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