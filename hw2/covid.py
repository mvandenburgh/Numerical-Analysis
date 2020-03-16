import numpy as np
import matplotlib.pyplot as plt
from numpy import sin, cos, sqrt, log, pi
from math import fmod
from random import random

x0 = 0.375205
x1 = 0.376047
def rng():
    global x0
    global x1

    temp = x0
    x0 = x1
    if x1 + temp > 1.0:
        x1 = fmod((x1 + temp) - 1.0, 1.0)    
    else:
        x1 = fmod((x1 + temp), 1.0)
    return x1


def box_muller_transform(mu, sigma):
    
    x1 = rng()
    x2 = rng()
    
    g1 = sqrt(-2 * log(x1)) * cos(2 * pi * x2)
    g2 = sqrt(-2 * log(x1)) * sin(2 * pi * x2)
    
    z1 = mu + g1 * sigma
    z2 = mu + g2 * sigma

    return z1, z2

dist = []

for i in range(0, 100000):
    dist.append(box_muller_transform(0.18, sqrt(0.08))[0])

plt.hist(dist, 1000, density=True)
plt.show()