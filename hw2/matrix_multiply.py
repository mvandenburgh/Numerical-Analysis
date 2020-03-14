# N = 300, 3.84 seconds
# N = 500, 14.69 seconds
# N = 1024, 163.27 seconds


from array import *
from random import uniform
from time import time

startTime = time()

N = 1024


def printMatrix(matrix):
    for i in matrix:
        for j in i:
            print(j, end = " ")
        print("")


def genMatrix():
    matrix = [[]]
    for i in range(0, N):
        row = []
        for j in range(0, N):
            row.insert(j, uniform(0, 2) - 1)
        matrix.insert(i, row)
    return matrix

def naiveMatrixMultiply(A,B):
    C = [[]]
    for i in range(0, N):
        row = []
        for j in range(0, N):
            sum = 0
            for k in range(0,N):
                sum = sum + A[i][k] * B[k][j]
            row.append(sum)
        C.append(row)
    return C


A = genMatrix()
B = genMatrix()
C = naiveMatrixMultiply(A,B)
runningTime = time() - startTime
printMatrix(C)

print()
print("Running time of program:", runningTime)


# double randomNumber() {
#     double random_value;

#     random_value = ((double)rand()/(double)(RAND_MAX/2)) - 1.0;

#     return random_value;
# }

# void naiveMatrixMultiply(double A[N][N], double B[N][N], double C[N][N]) {
#     for (int i = 0; i < N; i++) {
#         for (int j = 0; j < N; j++) {
#             double sum = 0;
#             for (int k = 0; k < N; k++) {
#                 sum = sum + A[i][k] * B[k][j];
#             }
#             C[i][j] = sum;
#         }
#     }
# }