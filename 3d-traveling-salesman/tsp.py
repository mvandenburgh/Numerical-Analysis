from math import fmod, sqrt, log, cos, sin, pi, factorial
import matplotlib.pyplot as plt
import itertools
from time import time

seeds = [0.123237445, ((time() % 60) % 1)]

def box_muller_transform(mu, sigma):
    x1 = uniform_random(0, 1)
    x2 = uniform_random(0, 1)

    g1 = sqrt(-2 * log(x1)) * cos(2 * pi * x2)
    g2 = sqrt(-2 * log(x1)) * sin(2 * pi * x2)
    
    z1 = mu + g1 * sigma
    z2 = mu + g2 * sigma
    return z1, z2

# Algorithm source: Lecture 2 notes
def uniform_random(a, b):
    temp = seeds[0]
    seeds[0] = seeds[1]
    if seeds[1] + temp > 1.0:
        seeds[1] = fmod((seeds[1] + temp) - 1.0, 1.0)  
    else:
        seeds[1] = fmod((seeds[1] + temp), 1.0)
    return a + (b-a) * seeds[1]

# computes the distance between 2 three dimensional points
def distance(p1, p2):
    # return ((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)**0.5
    return ((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2 + (p2[2] - p1[2])**2)**0.5


N = 9 # number of trees

x = [uniform_random(0, 100) for i in range(N)]
y = [uniform_random(0, 100) for i in range(N)]
z = [box_muller_transform(10, 1.777)[0] for i in range(N)]
coordinates = [(x,y,z) for (x,y,z) in zip(x,y,z)]


max_coord = 0 # index of maximum coordinate
max_coord_val = 0 # x1 + y1 + z1 of max coord
for i in range(N):
    if sum(list(coordinates[i])) > max_coord_val:
        max_coord_val = sum(list(coordinates[i]))
        max_coord = i

permutations = list(itertools.permutations(coordinates))

shortest_path_length = -1
shortest_path = [(-1,-1,-1) for i in range(N)]

for i in range(factorial(N)):
    if (permutations[i][0] != coordinates[max_coord]):
        continue # skip if it doesn't start with max_coord

    current_path_length = 0
    for j in range(len(permutations[i])-1):
        current_path_length += distance(permutations[i][j], permutations[i][j+1])

    if current_path_length < shortest_path_length or shortest_path_length == -1:
        shortest_path_length = current_path_length
        shortest_path = permutations[i]

print(coordinates)
print()
print(shortest_path)
print("Length:", shortest_path_length)