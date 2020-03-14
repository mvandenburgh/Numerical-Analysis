import numpy as np
import matplotlib.pyplot as plt
from random import random
from numpy import sin, cos, sqrt, log, pi

def genNumbers(count):
    vals = []
    for i in range(0, count):
        vals.append(random())
    return vals


def box_muller_transform(mu, sigma):
    
    x1 = random()
    x2 = random()
    
    g1 = sqrt(-2 * log(x1)) * cos(2 * pi * x2)
    g2 = sqrt(-2 * log(x1)) * sin(2 * pi * x2)
    
    z1 = mu + g1 * sigma
    z2 = mu + g2 * sigma

    return z1, z2

dist = []

for i in range(0, 1000000):
    dist.append(box_muller_transform(0.18, sqrt(0.08))[0])

plt.hist(dist, 100, density=True)
plt.show()