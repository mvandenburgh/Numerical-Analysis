from math import sqrt
import sys

def circleY(x, h, k, r):
    return (k + sqrt(r**2 - (x-h)**2), k - sqrt(r**2 - (x-h)**2))

def circleX(y, h, k, r):
    return (h + sqrt(r**2 - (y-k)**2), h - sqrt(r**2 - (y-k)**2))

def quadratic(a, b, c):
    d = b**2-4*a*c 
    # print(d)
    return (((-b + sqrt(d)) / 2 * a), ((-b - sqrt(d)) / 2 * a))

def heartEq(x):
    return quadratic(1, -(2*sqrt(abs(x))), -(-abs(x)-(x)**2+2))

# print(heartEq(0.5))
# sys.exit()


# center = - sqrt(2)
y = 1.414213562373095
keepGoing = True
center = y-0.0000001
while(center > -sqrt(2) and keepGoing):
    radius = y - center
    print("Center @ (0, " + str(center) + ")")
    print("Radius: " + str(radius))
    # current = -1.41421356 # left most point on heart
    current = center - radius
    while (current < 0):
        print("X Point: " + str(current) + "\n")
        if (heartEq(current)[1] > circleY(current, 0, center, radius)[1]):
            
            print("FINAL radius: " + str(radius))
            print("Final center: (0, " + str(center) + ")")
            keepGoing = False
            break
        current += 0.01
    center -= 0.0001
    # print(center)
    # while (current > center - radius): # while we're stil inside circle

# print(circleX(0, 0, 0, sqrt(2)))


x = -1.41421356
# while(True):
#     print('X = ' + str(x))
#     print(quadratic(1, -(2*sqrt(abs(x))), -(-abs(x)-(x)**2+2)))
#     # print('X = ' + str(x))
#     x-=0.000000001
# while (True)
# print(quadratic(1, -sqrt(2), -1.25))
print()
# print(quadratic(1, (2 * sqrt(0.5)), (-0.5 - (0.5**2) + 2)))

    
