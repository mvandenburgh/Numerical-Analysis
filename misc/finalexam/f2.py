import matplotlib.pyplot as plt
from numpy import sqrt, exp, log

def RK4(dydx, x0, y0, h, xn):
    """Classical 4th order Runge-Kutta method to solve first order differential equations numerically"""
    xs = [x0]
    ys = [y0]

    x = x0
    y = y0
    if x0 > xn: h *= -1

    while (True):
        k_0 = dydx(x, y)
        k_1 = dydx(x + h/2, y + h/2 * k_0)
        k_2 = dydx(x + h/2, y + h/2 * k_1)
        k_3 = dydx(x + h, y + h * k_2)
        k = (1/6) * (k_0 + 2.0*k_1 + 2.0*k_2 + k_3)

        x = x + h
        y = y + h * k
        if (x0 > xn and x < xn):
            break
        elif (x0 < xn and x > xn):
            break
        xs.append(x)
        ys.append(y)
    
    return xs, ys

dydx = lambda x, y: log(x**3 + y**3)

xs, ys = RK4(dydx, 0, 1, 0.01, 1)

plt.plot(xs, ys)
plt.show()