# AMS 326 Test 2
# 3/31/2020
# Michael VanDenburgh

import matplotlib.pyplot as plt
from numpy import sqrt, append, linalg

# Function to perform the classical Runge-Kutta method.
# http://people.bu.edu/andasari/courses/numericalpython/Week5Lecture9/PythonFiles/ex3_RK4thOrder_Numpy.py
# was used as a reference and some code/ideas were taken from there.
def RK4(dydx, yinit, x_range, h):
    x = x_range[0]
    y = yinit

    xsol = [x]
    ysol = [y]

    iterations = round((x_range[-1] - x_range[0])/h)
    for i in range(iterations):
        k1 = dydx(x ,y)

        y2= y + [i * h/2 for i in k1]

        k2 = dydx(x+h/2, y2)

        y3 = y + [i * h/2 for i in k2]

        k3 = dydx(x+h/2, y3)

        y4 = y + [i * h for i in k3]

        k4 = dydx(x+h, y4)

        for j in range(len(yinit)): # calculate slope
            y[j] = y[j] + (h/6)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j])

        x = x + h
        xsol = append(xsol, x)

        for i in range(len(y)):
            ysol = append(ysol, y[i])

    return [xsol, ysol]


# Represents y' function for use in RK4.
def dydx(x, y):
    dy = [float(0) for i in range(len(y))] # initialize an array of zeroes
    dy[0] = x + y[0] + x * y[0]
    return dy


# Calculate coefficients for interpolating polynomial
def get_polynomial_coefficients(xvals, yvals):
    A = [[x**i for i in range(len(xvals) -1, -1, -1)] for x in xvals]
    B = [y for y in yvals]
    coefficients = linalg.inv(A).dot(B) # solve system
    return coefficients


# Calculates points for interpolation graph based on
# provided coefficients
def get_interpolated_points(coefficients, xmin, xmax):
    degree = len(coefficients) - 1
    x = range(xmin, xmax)
    y = []

    for t in range(xmin, xmax):
        current = 0
        for a, i in zip(coefficients, range(degree, -1, -1)):
            current += (a * t**i)
        y.append(current)
    
    return x, y


# Returns index that value x is located at in array arr
def indexOf(x, arr):
    for i in range(0, len(arr)):
        if x == round(arr[i],2): # round to account for slight floating point difference
            return i
    return -1


# Part 1
x, y = RK4(dydx, [1], [0, 0.5], 0.01)

for i, j in zip(x, y):
    print(round(i,2), j) # print x/y values

# Part 2
plt.plot(x, y, 'r')
plt.xlabel('x')
plt.ylabel('y')
plt.show()


# Part 3

xvals = [0.1, 0.2, 0.3, 0.4, 0.5]
yvals = [y[indexOf(i, x)] for i in xvals]

poly_coeffcients = get_polynomial_coefficients(xvals, yvals)

print(poly_coeffcients)