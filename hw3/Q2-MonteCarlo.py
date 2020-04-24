import numpy as np
import matplotlib.pyplot as plt

from random import random

w = 1 # distance between parallel lines

d = 0.75 # diameter of disc

#samples = 1984444444 # number of discs to drop
samples = 1000

def toss_disk():
    distance = np.random.uniform(0, w/2) # distance from center of disc to nearest line
    
    if distance < d/2: # if the disc crosses a parallel line
        return 1
    
    else:
        return 0

observations = np.zeros(samples)

for i in range(samples):
    observations[i] = toss_disk()

# observations = [toss_disk() for i in range(samples)]


probabilities = [np.average(observations[0:i]) for i in range(samples)]

# plt.plot(range(samples), probabilities)

# plt.show()

print(probabilities)