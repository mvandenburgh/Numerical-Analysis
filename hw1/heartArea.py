##################################################
# This program calculates the area of the heart. #
##################################################

from math import sqrt
from time import time
startTime = time()


def integrate(f, a, b, n=10000000):
    """ Computes the definite integral of the function f from a to b using Simpson's Rule. 
    Optional argument n to provide number of subdivisions (default 10000).
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


def topLeft(x):
    """ Equation of the top part of the heart on the left side of the y-axis. """
    return sqrt(-x) + sqrt(-(x**2) - x - abs(x) + 2)

def bottomLeft(x):
    """ Equation of the bottom part of the heart on the left side of the y-axis. """
    return sqrt(-x) - sqrt(-(x**2) - x - abs(x) + 2)


area = 2 * (integrate(topLeft, -sqrt(2), 0) - integrate(bottomLeft, -sqrt(2), -1) + integrate(bottomLeft, 0, -1))
print("Area of heart = " + str(area))

print("Program finished in " + str(time()-startTime) + " seconds")