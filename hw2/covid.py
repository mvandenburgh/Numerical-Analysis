import numpy as np
import matplotlib.pyplot as plt
from numpy import sin, cos, sqrt, log, pi
from math import fmod
from random import random

x0 = 0.775205
x1 = 0.42942
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
    # x1 = rng()
    # x2 = rng()

    x1 = random()
    x2 = random()

    g1 = sqrt(-2 * log(x1)) * cos(2 * pi * x2)
    g2 = sqrt(-2 * log(x1)) * sin(2 * pi * x2)
    
    z1 = mu + g1 * sigma
    z2 = mu + g2 * sigma
    return z1, z2


x = []
y = [10]

for i in range(0, 45):
    x.append(box_muller_transform(0.18, sqrt(0.08))[0])

for i in range(0, 44):
    y.append(y[i] * (1 + x[i]))

plt.hist(x, 100, density=True)
plt.plot(range(0, 45), y)
plt.show()


for i in range(0, 45):
    x.append(box_muller_transform(-0.24, sqrt(0.04))[0])

for i in range(0, 45):
    y.append(y[i+44] * (1 + x[i+44]))

plt.plot(range(0, 90), y)
plt.show()