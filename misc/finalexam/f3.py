from time import time
from numpy import fmod, sqrt, log, cos, sin, pi, exp
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

seed0 = 0.775205 # first seed for pseudo-RNG
seed1 = ((time() % 60) % 1) # base second seed on current millisecond


def rng(a,b):
    """ Generate a random number drawn from a uniform distribution between a and b """
    global seed0, seed1
    temp = seed0
    seed0 = seed1
    if seed1 + temp > 1.0:
        seed1 = fmod((seed1 + temp) - 1.0, 1.0)  
    else:
        seed1 = fmod((seed1 + temp), 1.0)
    return a + (seed1 * (b-a))

def box_muller_transform(mu, sigma):
    """ Generate a random number drawn from a normal distribution with mean mu and variance sigma. """
    x1 = rng(0,1)
    x2 = rng(0,1)

    g1 = sqrt(-2 * log(x1)) * cos(2 * pi * x2)
    g2 = sqrt(-2 * log(x1)) * sin(2 * pi * x2)

    z1 = mu + g1 * sigma
    z2 = mu + g2 * sigma
    return z1

n1 = 60
n2 = 190
n = n1 + n2

mu1 = 0.005
sigma1 = 0.015
mu2 = 0.010
sigma2 = 0.013

index_daily_change_rate = [24600] # set Dow(0)

# Compute following 60 trading days
for i in range(n1-1):
    index_daily_change_rate.append(index_daily_change_rate[i] * (1 + box_muller_transform(mu1, sigma1)))

for i in range(n2):
    index_daily_change_rate.append(index_daily_change_rate[i] * (1 + box_muller_transform(mu2, sigma2)))

plt.plot(range(n), index_daily_change_rate)
plt.xlabel("Day")
plt.ylabel("Index change")
plt.title("Dow index daily values for next 60 + 190 days")
plt.show()

dow = lambda t, d0, beta: d0 * exp(beta * t)

params = curve_fit(f=dow, xdata=range(n1)[0:60:10], ydata=index_daily_change_rate[0:60:10])[0]

d0 = params[0]
beta = params[1]

print(params)

# ys = [dow(x, d0, beta) for x in range(n1)]
# plt.plot(range(n1), ys, color='r')
# plt.yscale = "logit"
# plt.show()