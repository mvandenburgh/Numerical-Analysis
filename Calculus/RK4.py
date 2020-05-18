"""Four_Step_Runge_Kutta_ODE1.py 
Implementation of the classic fourth-order method also refered as the
"original" Runge–Kutta method. This method is an implicit four step
Runge-Kutta method which solves an intial value problem numerically. 
"""

import matplotlib.pyplot as plt
from math import sqrt

def RK4(dydx, x0, y0, h, xn):
    xs = [x0]
    ys = [y0]

    x = x0
    y = y0
    if x0 > xn: h *= -1

    while (True):
        if (x0 > xn and x < xn):
            break
        elif (x0 < xn and x > xn):
            break
        k_0 = dydx(x, y)
        k_1 = dydx(x + h/2, y + h/2 * k_0)
        k_2 = dydx(x + h/2, y + h/2 * k_1)
        k_3 = dydx(x + h, y + h * k_2)
        k = (1/6) * (k_0 + 2.0*k_1 + 2.0*k_2 + k_3)

        x = x + h
        y = y + h * k
        xs.append(x)
        ys.append(y)
    
    return xs, ys

# Define w(x) and dy/dx to set up initial value problem
w = lambda x: (4*v0) * ( (x/a) - ( (x**2)/(a**2) ) )

# The function dy/dx. Error handling implemented if the function is undefined.
dydx = lambda x, y: ((y/x) - (w(x) / (vB * (x / (sqrt(x**2 + y**2)))) ) 
                    if (vB * (x / (sqrt(x**2 + y**2)))) != 0 else 0)

v0 = 14 # initial velocity
a = 7777 # river width
yA = 0 # y(a) = 0

h = 0.01
xn = 0

vB = 7
sol = RK4(dydx, 7777, 1, h, xn)

plt.plot(sol[0],sol[1])

vB = 14
sol = RK4(dydx, 7777, 1, h, xn)
plt.plot(sol[0],sol[1])

vB = 21
sol = RK4(dydx, 7777, 1, h, xn)
plt.plot(sol[0],sol[1])

plt.show()