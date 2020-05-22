import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

# 5.3 - 3D Traveling Salesman Problem

def distance(x1, y1, z1, x2, y2, z2):
    return ((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)**0.5

N = 5

x = [np.random.uniform(0, 10) for i in range(N)]
y = [np.random.uniform(0, 20) for i in range(N)]
z = [np.random.uniform(0, 30) for i in range(N)]

used = [0 for i in range(N)]

# ax.scatter(x,y,z,marker='o')

at = 0
point = at+1

for i in range(N):
    # while(used[at] == 1):
    #     at += 1
    #     if (at == N):
    #         at = 0

    while(point != at and used[point] != 1):
        point += 1
        if (point == N):
            point = 0
        
    current = point+1
    for j in range(N-1):
        if point >= N:
            point = 0
        if current >= N:
            current = 0
        if used[current] == 1 or current == at:
            continue
            
        if distance(x[point], y[point], z[point], x[at], y[at], z[at]) > distance(x[current], y[current], z[current], x[at], y[at], z[at]) and current != at:
            point = current
            
        # print(distance(x[point], y[point], z[point], x[at], y[at], z[at]))
        # print(distance(x[current], y[current], z[current], x[at], y[at], z[at]))
        current += 1
    print(at)
    print(point)
    print("\n")
    used[at] = 1
    at = point
    
    plt.plot([x[at], x[point]], [y[at],y[point]], [z[at],z[point]])

# for i in range(N):
#     print(x[i],y[i],z[i])

plt.show()