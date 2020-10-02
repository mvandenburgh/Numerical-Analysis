# Finds roots of function f given interval a and b using the bisection method
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

# Finds roots of function f given fPrime (its derivative)
# and an initial guess x0 using Newton's Method.
def newton(f, fPrime, x0, iterations=100000):
    for _ in range(iterations):
        x0 = x0 - f(x0) / fPrime(x0)
    return x0
