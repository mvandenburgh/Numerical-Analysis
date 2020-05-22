from datetime import datetime
from math import sqrt, log, cos, sin, pi, fmod
import matplotlib.pyplot as plt

seeds = [0, 0] # seeds for pseudo RNG

def init_seeds():
    seeds[0] = datetime.now().timestamp()
    while(seeds[0] > 1):
        seeds[0] /= 10
    seeds[1] = float(str(datetime.now().year) + str(datetime.now().month) + str(datetime.now().hour) + str(datetime.now().microsecond))
    while(seeds[1] > 1):
        seeds[1] /= 10

def uniform_random(a, b):
    """ Generate uniformally distributed pseudo random numbers
    between a and b. init_seeds() must be called before using.
    Parameters:
        a: lower limit
        b: upper limit
    """
    temp = seeds[0]
    seeds[0] = seeds[1]
    if seeds[1] + temp > 1.0:
        seeds[1] = fmod((seeds[1] + temp) - 1.0, 1.0)  
    else:
        seeds[1] = fmod((seeds[1] + temp), 1.0)
    return a + (b-a) * seeds[1]


def box_muller_transform(mu, sigma):
    """ Generate normally distributed random numbers from uniform 
    random numbers using the Box-Muller transform.
    Parameters:
        mu: mean of desired normal distribution
        sigma: square root of the variance of desired normal distribution
    """
    x1 = uniform_random(0, 1)
    x2 = uniform_random(0, 1)

    g1 = sqrt(-2 * log(x1)) * cos(2 * pi * x2)
    g2 = sqrt(-2 * log(x1)) * sin(2 * pi * x2)
    
    z1 = mu + g1 * sigma
    z2 = mu + g2 * sigma
    return z1

# Example code to generate normal distribution with mu=0 and sigma=1

init_seeds()

n = 1000000
normal_distribution = [box_muller_transform(0, 1) for i in range(n)]

plt.hist(normal_distribution, bins= int(n / 10000))
plt.show()