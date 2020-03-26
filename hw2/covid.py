import matplotlib.pyplot as plt
from numpy import sin, cos, sqrt, log, pi
from scipy.interpolate import lagrange
from scipy.optimize import curve_fit
from math import fmod, ceil
from random import random
from time import time

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

    # x1 = random()
    # x2 = random()

    g1 = sqrt(-2 * log(x1)) * cos(2 * pi * x2)
    g2 = sqrt(-2 * log(x1)) * sin(2 * pi * x2)
    
    z1 = mu + g1 * sigma
    z2 = mu + g2 * sigma
    return z1, z2


x = []
y = [10] # 10 infected patients on day 1

# generate normal distribution for first 45 days with mu=0.18 and sigma=0.08
for i in range(0, 45):
    x.append(box_muller_transform(0.18, sqrt(0.08))[0])


# calculate number of infected patients on each day
for i in range(0, 44):
    y.append(ceil(y[i] * (1 + x[i]))) # round up number of infected people to an integer


# plot number of infected patients in first 45 days
plt.plot(range(0, 45), y) 
plt.title(label="COVID-19 Infections for First 45 Days")
plt.show()


# generate normal distribution for second 45 days with mu=-0.24 and sigma=0.04
for i in range(0, 45):
    x.append(box_muller_transform(-0.24, sqrt(0.04))[0])


# calculate number of infected patients on each day
for i in range(0, 45):
    y.append(ceil(y[i+44] * (1 + x[i+44])))


# plot number of infected patients in second 45 days
plt.plot(range(0, 45), y[45:])
plt.title(label="COVID-19 Infections for Second 45 Days")
plt.show()

plt.plot(range(0, 90), y)
plt.title(label="COVID-19 Infections for First 90 Days")
plt.show()

poly = lagrange([x[9], x[18], x[27], x[36], x[45]], [y[9],y[18],y[27],y[36],y[45]]).c
# poly = lagrange([x[1], x[9], x[18], x[27], x[36], x[45]], [y[1], y[9], y[18], y[27], y[36], y[45]]).c

y_interpolated = []

for i in range(1, 46):
    current_point = 0
    for j in range(0, len(poly)):
        current_point += (i**(len(poly)-j)) * poly[j]
    y_interpolated.append(current_point)


plt.plot(range(0,45), y_interpolated)
plt.title("COVID-19 Infections for First 45 Days (Interpolated from " + str(len(poly)) + " points)")
plt.show()

# curve_fit(tck, xdata=x_points, ydata=y_points)

def show_uniform_distribution(n):
    dist = []
    for i in range(0,n):
        dist.append(rng())
    plt.hist(dist)
    plt.show()