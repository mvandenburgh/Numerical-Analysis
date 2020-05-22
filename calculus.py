def integrate(f, a, b, n=1000000):
    """ Approximates the definite integral of the function f from a to b using Simpson's Rule. 
        Exception handling present to handle slight round-off/truncation errors that cause the function to be undefined.
        Parameters:
            f: function to integrate
            a: lower limit of integration
            b: upper limit of integration
            n: max number of subdivisions
    """
    deltaX = (b-a) / n
    sum = 0
    for i in range(n):
        x = a + i * deltaX
        y = f(x)
        if i == 0:
            try:
                sum += y
            except:
                pass
        elif i % 2 == 1:
            try: 
                sum += (4 * y)
            except:
                pass
        else:
            try: 
                sum += (2 * y)
            except:
                pass
    sum += f(a + n * deltaX)
    return (deltaX / 3) * sum


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