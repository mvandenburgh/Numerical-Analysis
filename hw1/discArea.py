from math import sqrt
from time import time
startTime = time()

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

incrementBy = 0.0001 # value to increment by when checking for intersection of disc and heart.
y=sqrt(2)
keepGoing = True # boolean flag to terminate while loop when radius/center are found
# center = 0.3824 Use this to get closer
center = y - incrementBy
while(center > -sqrt(2) and keepGoing): # while the center of the disc is inside the heart
    radius = y - center
    current = 0 - radius
    while (current < 0):
        if (heartEq(current + incrementBy)[1] >= discEq(current + incrementBy, 0, center, radius)[1]):
            print("Approximated radius disc: " + str(radius))
            print("Approximated center of disc (x, y): (0, " + str(center) + ")")
            print("Program finished in " + str(time() - startTime) + " seconds.")
            keepGoing = False
            break
        current += incrementBy
    center -= incrementBy
