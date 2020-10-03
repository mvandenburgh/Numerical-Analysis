from calculus import differentiate

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


# Finds roots of function f, given an initial guess x0, using Newton's Method.
def newton(f, x0, iterations=100000):
    for _ in range(iterations):
        x0 = x0 - f(x0) / differentiate(f, x0)
    return x0


# Finds roots of function f given two initial guesses,
# x0 and x1, using the secant method
def secant(f, x0, x1, iterations=100000):
    for _ in range(iterations):
        if (f(x1) - f(x0)) == 0:
            break
        temp = x1
        x1 = x1 - (f(x1) * (x1 - x0)) / (f(x1) - f(x0))
        x0 = temp
    return x1
