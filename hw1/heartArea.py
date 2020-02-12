from math import sqrt
from scipy import integrate
from time import time
startTime = time()


# equation of the top part of the heart on the left side of the y-axis
def topLeft(x):
    return sqrt(-x) + sqrt(-(x**2) - x - abs(x) + 2)

# equation of the bottom part of the heart on the left side of the y-axis
def bottomLeft(x):
    return sqrt(-x) - sqrt(-(x**2) - x - abs(x) + 2)

area = 2 * (integrate.quad(topLeft, -sqrt(2), 0)[0] - integrate.quad(bottomLeft, -sqrt(2), -1)[0] + integrate.quad(bottomLeft, 0, -1)[0])
print("Area of heart: " + str(area))

print("Program finished in " + str(time()-startTime) + " seconds")