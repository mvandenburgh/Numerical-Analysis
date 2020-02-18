from math import cos

# Finds roots of functino f given interval a and b using the bisection method
def bisection(f, a, b, epsilon=0.00000000001, n=10000000000):
    for i in range(n):
        x = (a + b) / 2

        if f(x) == 0 or (b - a) / 2 < epsilon:
            return x

        if f(x) * f(a) > 0:
            a = x
        elif f(x) * f(b) > 0:
            b = x
        else:
            print("Bisection method failed...")
            return None

