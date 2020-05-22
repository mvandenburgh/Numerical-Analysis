from time import time
from numpy import fmod, sqrt, log, cos, sin, pi
import numpy as np
from pylab import show,hist,subplot,figure
import matplotlib.pyplot as plt

seed0 = 0.775205 # first seed for pseudo-RNG
seed1 = ((time() % 60) % 1) # base second seed on current millisecond
def rng():
    global seed0
    global seed1

    temp = seed0
    seed0 = seed1
    if seed1 + temp > 1.0:
        seed1 = fmod((seed1 + temp) - 1.0, 1.0)  
    else:
        seed1 = fmod((seed1 + temp), 1.0)
    return seed1


def box_muller_transform(mu, sigma):
    x1 = rng()
    x2 = rng()

    g1 = sqrt(-2 * log(x1)) * cos(2 * pi * x2)
    g2 = sqrt(-2 * log(x1)) * sin(2 * pi * x2)
    
    z1 = mu + g1 * sigma
    z2 = mu + g2 * sigma
    return z1


n = 12
matrix = [[box_muller_transform(0,1) for i in range(n)] for i in range(n)]

vals, vects = np.linalg.eig(matrix)
maxcol = list(vals).index(max(vals))
eigenvect = vects[:,maxcol]
eigenval = max(vals)
print(eigenvect)
print(eigenval)


bins = []
for i in range(n):
    for j in range(n):
        bins.append(matrix[i][j])
hist(bins)
plt.show()