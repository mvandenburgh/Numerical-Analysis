import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Boat Trajectories")

# Define w(x) and dy/dx to set up initial value problem
w = lambda x: (4*v0) * ( (x/a) - ( (x**2)/(a**2) ) )

# The function dy/dx. Error handling implemented if the function is undefined.
dydx = lambda x, y: ((y/x) - (w(x) / (vB * (x / (np.sqrt(x**2 + y**2)))) ) 
                    if (vB * (x / (np.sqrt(x**2 + y**2)))) != 0 else 0)

v0 = 14 # initial velocity
a = 7777 # river width
yA = 0 # y(a) = 0


vB = 7
sol = solve_ivp(fun=dydx, t_span=[a, 0], y0=[0])
plt.plot(sol.t, sol.y[0], label="vB = " + str(vB))

vB = 14
sol = solve_ivp(fun=dydx, t_span=[a, 0], y0=[0])
plt.plot(sol.t, sol.y[0], label="vB = " + str(vB))

vB = 21
sol = solve_ivp(fun=dydx, t_span=[a, 0], y0=[0])
plt.plot(sol.t, sol.y[0], label="vB = " + str(vB))

plt.legend()
plt.show()
