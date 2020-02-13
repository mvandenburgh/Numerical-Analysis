from math import cos

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

f = lambda x: 2.020**(-x**3)-(x**3)*(cos(x**4))-1.984

print(bisection(f, -0.9, -0.8))
print(bisection(f, 1.2, 1.32))
print(bisection(f, 1.3, 1.5))
print(bisection(f, 1.65, 1.75))
print(bisection(f, 1.75, 1.85))
print(bisection(f, 1.9, 2.0))