from random import random

N = 3000000

bins = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

for i in range(0, N):
    x = random()
    if x > 0 and x < 0.1:
        bins[0] += 1
    elif x > 0.1 and x < 0.2:
        bins[1] += 1
    elif x > 0.2 and x < 0.3:
        bins[2] += 1
    elif x > 0.3 and x < 0.4:
        bins[3] += 1
    elif x > 0.4 and x < 0.5:
        bins[4] += 1
    elif x > 0.5 and x < 0.6:
        bins[5] += 1
    elif x > 0.6 and x < 0.7:
        bins[6] += 1
    elif x > 0.7 and x < 0.8:
        bins[7] += 1
    elif x > 0.8 and x < 0.9:
        bins[8] += 1
    else:
        bins[9] += 1

print("(0.0, 0.1):", bins[0])
print("(0.1, 0.2):", bins[1])
print("(0.2, 0.3):", bins[2])
print("(0.3, 0.4):", bins[3])
print("(0.4, 0.5):", bins[4])
print("(0.5, 0.6):", bins[5])
print("(0.6, 0.7):", bins[6])
print("(0.7, 0.8):", bins[7])
print("(0.8, 0.9):", bins[8])
print("(0.9, 1.0):", bins[9])


mu = 10
sigma = 5


from numpy import random, sqrt, log, sin, cos, pi
from pylab import show,hist,subplot,figure

# Box-Muller transformation function
def box_muller_transform(u1,u2, mu, sigma):
  z1 = mu + (sqrt(-2*log(u1))*cos(2*pi*u2)) * sigma
  z2 = mu + (sqrt(-2*log(u1))*sin(2*pi*u2)) * sigma
  return z1,z2


# Generate random values between 0 and 1
u1 = random.rand(N)
u2 = random.rand(N)

# run the Box-Muller transformation
z1,z2 = box_muller_transform(u1,u2, mu, sigma)

hist(z1)
show()