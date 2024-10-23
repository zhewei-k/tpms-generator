def f(x):
    return x**6 - x - 1     #FEA calculated elastic modulus

def df(x):
    return 6*x**5 - 1

#Implement Secant method

def Secant(x0, f, max_it, tol, count = None, x_prev = None):
    
    if count is None:
        count = 0
        x_vals = []
        x_prev = x0 + tol + 1
        
    if count > max_it:
        return "Failed to converge"
    

    grad = (f(x0)-f(x_prev))/(x0-x_prev)
    xn = x0 - f(x0)/grad
    x_vals.append(xn)

    if abs(x0 - x_prev) < tol:
        return x0
    else:
        return Secant(xn, f, max_it, tol, count + 1, x0)
    
def Newton(x0, f, df, max_it, tol, count = None, x_prev = None):
    if count is None:
        count = 0 
        x_prev = x0 + tol + 1 
        
    if count > max_it:
        return "Failed to converge"
    
    print(x0 - f(x0)/df(x0))
    
    if abs(x0 - x_prev) < tol:
        return x0
    
    else:
        return Newton(x0 - f(x0)/df(x0), f, df, max_it, tol, count + 1, x0)
    
print(Secant(1, f, 100, 1e-5))
print(Newton(1, f, df, 100, 1e-5))
