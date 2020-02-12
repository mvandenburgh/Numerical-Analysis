from math import sqrt

import numpy as np


def circleY(x, h, k, r):
    return (k + sqrt(r**2 - (x-h)**2), k - sqrt(r**2 - (x-h)**2))

def circleX(y, h, k, r):
    return (h + sqrt(r**2 - (y-k)**2), h - sqrt(r**2 - (y-k)**2))

def quadratic(a, b, c):
    d = b**2-4*a*c 
    # print(d)
    return (((-b + sqrt(d)) / 2 * a), ((-b - sqrt(d)) / 2 * a))

def heartEq(x):
    return quadratic(1, -(2*sqrt(abs(x))), -(-abs(x)-(x)**2+2))

def distance(xc, yc, xp, yp):
    return sqrt((xp-xc)**2 + (yp-yc)**2)

incrementBy = 0.00001
y=sqrt(2)
keepGoing = True
# center = 0.3824
center = y-0.0001

while(center > -sqrt(2) and keepGoing): # while the center of the disc is inside the heart
    radius = y - center
    print("Center @ (0, " + str(center) + ")")
    print("Radius: " + str(radius) + '\n')
    current = 0 - radius
    while (current < 0):
        if (heartEq(current + incrementBy)[1] >= circleY(current + incrementBy, 0, center, radius)[1]):
            print("FINAL radius: " + str(radius))
            print("Final center: (0, " + str(center) + ")")
            keepGoing = False
            break
        current += incrementBy
    center -= incrementBy
