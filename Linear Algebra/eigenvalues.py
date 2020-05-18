import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
import scipy.optimize

n = 5
a = [
    [3.0000, 6.0000, 7.0000, 5.0000, 3.0000],
    [5.0000, 6.0000, 2.0000, 9.0000, 1.0000],
    [2.0000, 7.0000, 0.0000, 9.0000, 3.0000],
    [6.0000, 0.0000, 6.0000, 2.0000, 6.0000],
    [1.0000, 8.0000, 7.0000, 9.0000, 2.0000]
]

# Finds roots of functino f given interval a and b using the bisection method
def bisection(f, a, b, epsilon=0.00000000001, n=10000000000):
    for i in range(n):
        x = (a + b) / 2
        y = f(x)
        if y == 0 or (b - a) / 2 < epsilon:
            return x

        if y * f(a) > 0:
            a = x
        elif y * f(b) > 0:
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
    print(lam)
    # mat = deepcopy(a)
    for i in range(n):
        a[i][i] -= lam
    ret = determinant(a)
    for i in range(n):
        a[i][i] += lam
    return ret

def graph_function(f, a, b, h=0.001):
    xs = []
    ys = []
    x = a
    while(x < b):
        print(x)
        xs.append(x)
        ys.append(f(x))
        x += h
    plt.plot(xs, ys)
    plt.show()


# graph_function(char_poly, -4, 3, 0.01)

print("root:", bisection(char_poly, -4, -1))
# print("root:",bisection(char_poly, -1, 1))
# print(bisection(char_poly, 4, 6))
# print(bisection(char_poly, 22, 26))