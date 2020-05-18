from time import time
from math import fmod, sqrt, log, sin, cos, pi
from copy import deepcopy

import matplotlib.pyplot as plt

seed0 = 0.775205 # first seed for pseudo-RNG
seed1 = ((time() % 60) % 1) # base second seed on current millisecond


def rng(a,b):
    """ Generate a random number drawn from a uniform distribution between a and b """
    global seed0, seed1
    temp = seed0
    seed0 = seed1
    if seed1 + temp > 1.0:
        seed1 = fmod((seed1 + temp) - 1.0, 1.0)  
    else:
        seed1 = fmod((seed1 + temp), 1.0)
    return a + (seed1 * (b-a))

def box_muller_transform(mu, sigma):
    """ Generate a random number drawn from a normal distribution with mean mu and variance sigma. """
    x1 = rng(0,1)
    x2 = rng(0,1)

    g1 = sqrt(-2 * log(x1)) * cos(2 * pi * x2)
    g2 = sqrt(-2 * log(x1)) * sin(2 * pi * x2)

    z1 = mu + g1 * sigma
    z2 = mu + g2 * sigma
    return z1

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

n = 32

matrix = [[box_muller_transform(0, 1) for i in range(n)] for i in range(n)]


def char_poly(lam):
    global matrix, n
    mat = deepcopy(matrix)
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


graph_function(char_poly, -10, -9, 0.1)