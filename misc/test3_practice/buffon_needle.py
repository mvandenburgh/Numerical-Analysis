from math import fmod, sin, cos, pi

seeds = [0.12345, 0.43212]

line_width = 1 # distance between lines
needle_length = 0.75 # needle length

def uniform_random(a, b):
    temp = seeds[0]
    seeds[0] = seeds[1]
    if seeds[1] + temp > 1.0:
        seeds[1] = fmod((seeds[1] + temp) - 1.0, 1.0)  
    else:
        seeds[1] = fmod((seeds[1] + temp), 1.0)
    return a + (b-a) * seeds[1]

def drop_needle():
    distance = uniform_random(0, line_width/2)
    theta = uniform_random(0, pi/2)
    
    if distance < needle_length/2 * sin(theta):
        return 1
    else:
        return 0

def toss_triangle(h, line_width):
    ''' h - height of triangle '''
    distance = uniform_random(0, line_width/2)
    theta = uniform_random(0, pi/3)

    if distance < (h/2) * cos(theta):
        return 1
    else:
        return 0

def toss_square(diagonal_length, line_width):
    ''' diagonal_length - length of square's diagonal '''
    distance = uniform_random(0, line_width/2)
    theta = uniform_random(0, pi/4)

    if diagonal_length/2 * cos(theta) > distance:
        return 1
    else:
        return 0

def toss_hexagon(a, line_width):
    ''' a - length of hexagon edge '''
    distance = uniform_random(0, line_width/2)
    theta = uniform_random(0, pi/6)

    if distance < a * cos(theta):
        return 1
    else:
        return 0


def toss_octagon(a, line_width):
    ''' a - length of octagon edge '''
    distance = uniform_random(0, line_width/2)
    theta = uniform_random(0, pi/8)

    if distance < a * cos(theta):
        return 1
    else:
        return 0

def toss_decagon(a, line_width):
    ''' a - length of decagon edge '''
    distance = uniform_random(0, line_width/2)
    theta = uniform_random(0, pi/10)

    if distance < a * cos(theta):
        return 1
    else:
        return 0

def toss_12gon(a, line_width):
    ''' a - length of decagon edge '''
    distance = uniform_random(0, line_width/2)
    theta = uniform_random(0, pi/12)

    if distance < a * cos(theta):
        return 1
    else:
        return 0

def toss_ngon(n, a, line_width):
    ''' a - length of an edge '''
    distance = uniform_random(0, line_width/2)
    theta = uniform_random(0, pi/n)

    if distance < a * cos(theta):
        return 1
    else:
        return 0

samples = 10000000

sum = 0
for i in range(samples):
    # sum += toss_square(1, 1)
    # sum += toss_triangle(1,1)
    sum += toss_ngon(3, .5, 1)

print(sum / samples)