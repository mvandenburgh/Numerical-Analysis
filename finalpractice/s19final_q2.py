from math import sqrt, sin, pi, pow
from scipy import random

def f(x, y):
    return sqrt( ( pow((pow(x, 2) + pow(y, 2)), 3) ) / (4 * pow(x, 2)) )


def r(theta):
    return sin(2 * theta)

def integrate(f, a, b, n=1000000):
    """ Computes the definite integral of the function f from a to b using Simpson's Rule. 
    Optional argument n to provide number of subdivisions (default 1000000).
    Exception handling present to handle slight round-off/truncation errors that cause the function to be undefined."""
    deltaX = (b-a) / n
    sum = 0
    for i in range(n):
        x = a + i * deltaX
        if i == 0:
            try:
                sum += f(x)
            except:
                pass
        elif i % 2 == 1:
            try: 
                sum += (4 * f(x))
            except:
                pass
        else:
            try: 
                sum += (2 * f(x))
            except:
                pass
    sum += f(a + n * deltaX)
    return (deltaX / 3) * sum

max_x = 0.76942088429773

unit_circle = lambda x: sqrt(1 - x**2)


# limits of integration
a = 0
b = 2 * pi

area = integrate(lambda theta: (1/2) * r(theta)**2, a, b)
print(area)

N = 1000000 # number of points to generate for MC simulation

sum = sum([(1/2) * (r(random.uniform(a, b)))**2 for i in range(N)])

approx_area = (b-a)/float(N) * sum

print(approx_area)