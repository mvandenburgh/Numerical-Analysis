import math
import random

def randnumgen(x, y):
    return (x + y) % 1.0

def convert(x1, x2):
    g1 = math.sqrt(-2 * math.log(x1, math.e)) * math.cos(2 * math.pi * x2)
    g2 = math.sqrt(-2 * math.log(x1, math.e)) * math.sin(2 * math.pi * x2)
    return g1, g2

def genlist(seed1, seed2,iterations, mu, sigma):
    s1 = seed1
    s2 = seed2
    l = []
    for i in range(iterations):
        temp = s2
        s2 = randnumgen(s1, s2)
        s1 = temp
        g1, g2 = convert(s1, s2)
        l.append(round(mu + g1 * sigma, 5))
        l.append(round(mu + g2 * sigma, 5))
    return l

def genPoints():
    z = genlist(0.78237423, 0.2374832, 5, 10, 1.777)
    p = []
    s1 = 0.83495734
    s2 = 0.23784224
    s3 = 0.73485343
    s4 = 0.78345844
    for i in range(9):
        temp1 = s2
        s2 = randnumgen(s1, s2)
        s1 = temp1
        temp2 = s4
        s4 = randnumgen(s3, s4)
        s3 = temp2
        x = round(s2*100,5)
        y = round(s4*100,5)
        p.append([x,y,z[i]])
    return p


def dist(v1, v2):
    return math.sqrt((v1[0]-v2[0])**2 + (v1[1] - v2[1])**2 + (v1[2] - v2[2])**2)

def pathDist(points):
    total = 0
    for i in range(len(points) - 1):
        total += dist(points[i],points[i+1])
    return total

def pathDist2(points, a, b):
    total = 0
    for i in range(len(points) - 1):
        if i == a:
            index = b
        elif i == b:
            index = a
        else:
            index = i
        if (i+1) == a:
            index2 = b
        elif (i+1) == b:
            index2 = a
        else:
            index2 = i + 1
        total += dist(points[index], points[index2])
    return total

def toSwap():
    a = int(random.random()*8) + 1
    b = int(random.random()*8) + 1
    while a == b:
        b = int(random.random() * 8) + 1
    return a, b

def prob(ey, ex,t0, k, type):
    t = type(t0, k)
    try:
        p = math.e**(-(ey-ex)/t)
        return min(1,p)
    except Exception:
        return 0

def boltAnnealing(t0, k):
    try:
        return t0 / math.log(k)
    except ZeroDivisionError:
        return 0

def simulatedAnnealing(type, points):
    d = pathDist(points)
    counter = 0
    iteration = 1
    while counter < 100:
        a, b = toSwap()
        temp = pathDist2(points, a, b)
        p = prob(temp, d, 100, iteration, type)
        print("Current path distance:", d)
        chance = random.random()
        if p > chance:
            temp = points[a]
            points[a] = points[b]
            points[b] = temp
            d = pathDist(points)
            counter = 0
        else:
            counter += 1
        iteration += 1

def starting(points):
    max = 0
    index = 0
    for i in range(len(points)):
        p = points[i]
        temp = p[0] + p[1] + p[2]
        if temp > max:
            max = temp
            index = i
    return index

def divide(list):
    x = []
    y = []
    z = []
    for point in list:
        x.append(point[0])
        y.append(point[1])
        z.append(point[2])
    return x, y, z

if __name__ == "__main__":
    # Seeds can be changed at genPoints()
    # Generate 9 points
    points = genPoints()

    # Put max x+y+z first point
    index = starting(points)
    temp = points[index]
    points[index] = points[0]
    points[0] = temp

    # Starting Path
    print("Starting Path", points)
    print("Path total distance:", pathDist(points))

    # Simulated Annealing
    simulatedAnnealing(boltAnnealing, points)

    # Final Path
    print("Optimal-distance path", points)
    print("Path total distance:", pathDist(points))


