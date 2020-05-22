import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos, sqrt, log, pi, exp
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


# graph number of infected patients in first 45 days
plt.scatter(range(0, 45), y) 
plt.title(label="COVID-19 Infections for First 45 Days")
plt.xlabel("Days")
plt.ylabel("Infected patients")
plt.show()


# generate normal distribution for second 45 days with mu=-0.24 and sigma=0.04
for i in range(0, 45):
    x.append(box_muller_transform(-0.24, sqrt(0.04))[0])


# calculate number of infected patients on each day
for i in range(0, 45):
    y.append(ceil(y[i+44] * (1 + x[i+44])))


# graph number of infected patients in second 45 days
plt.scatter(range(45, 90), y[45:])
plt.title(label="COVID-19 Infections for Second 45 Days")
plt.xlabel("Days")
plt.ylabel("Infected patients")
plt.show()


# calculate coefficients for polynomial
A = [
    [9**4, 9**3, 9**2, 9**1, 9**0],
    [18**4, 18**3, 18**2, 18**1, 18**0],
    [27**4, 27**3, 27**2, 27**1, 27**0],
    [36**4, 36**3, 36**2, 36**1, 36**0],
    [45**4, 45**3, 45**2, 45**1, 45**0]
]

B = [y[9], y[18], y[27], y[36], y[45]]

poly = np.linalg.inv(A).dot(B) # coefficients of polynomial

# calculate y values for interpolation
y_interpolated = []
for i in range(1, 46):
    current_point = i**4 * poly[0] + i**3 * poly[1] + i**2 * poly[2] + i**1 * poly[3] + i**0 * poly[4]
    y_interpolated.append(current_point)



# plot interpolated curve and original
plt.plot(range(0, 45), y_interpolated, label="interpolated", color='r')
plt.scatter(range(0, 45), y[:45], label="original")
plt.title("COVID-19 Infections for First 45 Days (Interpolated from " + str(len(poly)) + " points)")
plt.xlabel("Days")
plt.ylabel("Infected patients")
plt.legend()
plt.show()


y_data_fitted = []

curve_fit_function = lambda t, a, b: a * exp(b * t)

params = curve_fit(f=curve_fit_function, xdata=x[45:], ydata=y[45:])[0]

# calculate y values for fitted curve
for i in range(0, 45):
    y_data_fitted.append(curve_fit_function(x[45:][i], params[0], params[1]))

plt.scatter(range(45, 90), y[45:], label="original", color='black')
plt.plot(range(45, 90), y_data_fitted, label="fitted", color='r')
plt.title("COVID-19 Infections for Second 45 Days (Curve Fitting)")
plt.xlabel("Days")
plt.yscale("symlog")
plt.ylabel("Infected patients")
plt.legend()
plt.show()


plt.title("All results")
plt.scatter(range(0, 45), y[:45], label="First 45 Days", color="r")
plt.scatter(range(45, 90), y[45:], label="Second 45 Days", color="g")
plt.plot(range(0, 45), y_interpolated, label="Interpolation", color='k')
plt.plot(range(45, 90), y_data_fitted, label="Curve Fitting", color="b")
plt.xlabel("Days")
plt.ylabel("Infected patients")
plt.legend()
plt.show()