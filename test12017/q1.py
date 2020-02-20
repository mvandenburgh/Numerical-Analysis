from math import cos, pi

def simpson(f, a, b, n=1000000):
    """ Computes the definite integral of the function f from a to b using Simpson's Rule. 
    Optional argument n to provide number of subdivisions (default 10000).
    Exception handling present to handle slight round-off/truncation errors that cause the function to be undefined."""
    deltaX = (b-a) / n
    sum = 0
    for i in range(n):
        x = a + i * deltaX
        if i == 0:
            try:
                sum += f(x)
            except:
                pass
        elif i % 2 == 1:
            try: 
                sum += (4 * f(x))
            except:
                pass
        else:
            try: 
                sum += (2 * f(x))
            except:
                pass
    sum += f(a + n * deltaX)
    return (deltaX / 3) * sum

f = lambda theta: (1 - cos(4*theta))

print((1/4) * simpson(f, 0, 2 * pi))
print(pi/2)