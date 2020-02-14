#################################################
# This program calculates the area of the disc. #
#################################################


from math import sqrt, pi
from time import time
startTime = time()

# Calculates the distance between two points
def distance(x1, y1, x2, y2):
    return sqrt((x2-x1)**2 + (y2-y1)**2)

# Quadratic equation solver
def quadratic(a, b, c):
    d = b**2-4*a*c
    return (((-b + sqrt(d)) / 2 * a), ((-b - sqrt(d)) / 2 * a))

# equation of the disc
def discEq(x, h, k, r):
    return (k + sqrt(r**2 - (x-h)**2), k - sqrt(r**2 - (x-h)**2))

# equation of the heart
def heartEq(x):
    return quadratic(1, -(2*sqrt(abs(x))), -(-abs(x)-(x)**2+2))

incrementBy = 0.001 # value to increment by when checking for intersection of disc and heart.
y=sqrt(2) # y-value of the top of the disc
keepGoing = True # boolean flag to terminate while loop when radius/center are found
center = 1 # approximate guess
iterations = 4 # increase to get more accurate results. Anything > 2 takes more than a second.
for i in range(iterations):
    keepGoing = True
    while(center > -sqrt(2) and keepGoing): # while the center of the disc is inside the heart
        radius = y - center
        current = 0 - radius # x value
        while (current < 0):
            if (distance(current, heartEq(current)[1], current, discEq(current, 0, center, radius)[1]) <= incrementBy):
                print('\nIteration #', i+1)
                print("Approximated radius disc: " + str(radius))
                print("Approximated center of disc (x, y): (0.0, " + str(center) + ")")
                keepGoing = False
                break
            current += incrementBy
        center -= incrementBy
    center += incrementBy * 10 # increase the approx center slightly to account for possible error
    incrementBy /= 10
    

print("\nArea of disc = " + str(pi * radius**2))
print("\nProgram finished in " + str(time() - startTime) + " seconds.")
