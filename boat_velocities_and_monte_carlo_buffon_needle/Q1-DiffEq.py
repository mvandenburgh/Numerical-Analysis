import matplotlib.pyplot as plt
from numpy import sqrt

def RK4(dydx, x0, y0, h, xn):
    """ Approximates a first order differential equation using the classical 4th order Runge-Kutta 
        Parameters:
            dydx: function representing differential equation to solve
            x0: initial x value
            y0: initial y value
            h: step-size
            xn: x-value stopping point
    """
    xs = [x0]
    ys = [y0]
    
    x = x0
    y = y0
    if x0 > xn: h *= -1

    while (True):
        k0 = dydx(x, y)
        k1 = dydx(x + h/2, y + h/2 * k0)
        k2 = dydx(x + h/2, y + h/2 * k1)
        k3 = dydx(x + h, y + h * k2)
        k = (1/6) * (k0 + 2.0*k1 + 2.0*k2 + k3)

        x = x + h
        y = y + h * k
        if (x0 > xn and x < xn):
            break
        elif (x0 < xn and x > xn):
            break
        xs.append(x)
        ys.append(y)
    
    return xs, ys

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Boat Trajectories")

# Define w(x) and dy/dx to set up initial value problem
w = lambda x: (4*v0) * ( (x/a) - ( (x**2)/(a**2) ) )

# The function dy/dx. Error handling implemented if the function is undefined.
dydx = lambda x, y: ((y/x) - (w(x) / (vB * (x / (sqrt(x**2 + y**2)))) ) 
                    if (vB * (x / (sqrt(x**2 + y**2)))) != 0 else 0)

v0 = 14 # initial velocity
a = 7777 # river width
yA = 0 # y(a) = 0

x0 = a # initial x value
y0 = 0 # initial y value
xn = 0 # stop at x=0

for vB in [7, 14, 21]:
    sol = RK4(dydx=dydx, x0=x0, y0=y0, xn=xn, h=0.05)
    plt.plot(sol[0], sol[1], label="vB = " + str(vB))

plt.legend()
plt.show()