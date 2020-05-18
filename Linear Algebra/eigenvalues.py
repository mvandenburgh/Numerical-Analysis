import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
import scipy.optimize

n = 5
a = [
    [1,2,3,4,5],
    [6,7,8,9,0],
    [1,2,3,4,5],
    [6,7,8,9,0],
    [9,8,7,6,5]
]
# \begin{pmatrix}1&2&3&4&5\\ 6&7&8&9&0\\ 1&2&3&4&5\\ 6&7&8&9&0\\ 9&8&7&6&5\end{pmatrix}

# Finds roots of functino f given interval a and b using the bisection method
def bisection(f, a, b, epsilon=0.00000000001, n=10000000000):
    for i in range(n):
        x = (a + b) / 2

        if f(x) == 0 or (b - a) / 2 < epsilon:
            return x

        if f(x) * f(a) > 0:
            a = x
        elif f(x) * f(b) > 0:
            b = x
        else:
            print("Bisection method failed...")
            return None

def determinant_helper(matrix, mul):
    size = len(matrix)
    if size == 1:
        return mul * matrix[0][0]
    else:
        sign = -1
        sum = 0
        for i in range(size):
            current = []
            for j in range(1, size):
                temp = []
                for k in range(size):
                    if k != i:
                        temp.append(matrix[j][k])
                current.append(temp)
            sign *= -1
            sum += mul * determinant_helper(current, sign * matrix[0][i])
        return sum

def determinant(matrix):
    return determinant_helper(matrix, 1)

def char_poly(lam):
    mat = deepcopy(a)
    for i in range(n):
        mat[i][i] -= lam
    return determinant(mat)

def graph_function(f, a, b, h=0.001):
    xs = []
    ys = []
    x = a
    while(x < b):
        xs.append(x)
        ys.append(f(x))
        x += h
    plt.plot(xs, ys)
    plt.show()


graph_function(char_poly, 0, 5)

print(bisection(char_poly, -4, -3))
print(bisection(char_poly, -1, 1))
print(bisection(char_poly, 4, 6))
print(bisection(char_poly, 22, 26))