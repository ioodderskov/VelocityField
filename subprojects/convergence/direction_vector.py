from __future__ import division
import sympy as sym
import scipy as sp
#x, y, z, t = sym.symbols('x y z t')
#k, m, n = sym.symbols('k m n', integer=True)
#f, g, h = sym.symbols('f g h', cls=Function)
#
sym.init_printing(use_unicode=True)
#print sym.diff(x**2/2, x)




a = sym.symbols('a', positive=True)
B = sym.symbols('B')
r, theta = sym.symbols('r theta')
vr, vtheta, vphi = sym.symbols('vr vtheta vphi')
ut, ur, utheta, uphi = sym.symbols('ut ur utheta uphi')
kt, kr = sym.symbols('kt kr')

alpha = sym.symbols('alpha')

# The components of the supernovas 4-velocity:
Ve = sym.sqrt(1-a**2*vr**2)

ut = 1/Ve
ur = vr/Ve
utheta = vtheta/Ve
uphi = vphi/Ve


alpha = (1/a + vr)/Ve
#
time_term = -(kt+alpha*ut)**2
radial_term = a**2 * (kr+alpha*ur)**2
theta_term = a**2 * r**2*(alpha*utheta)**2
phi_term = a**2 * r**2*sym.sin(theta)*(alpha*uphi)**2 
#
ndotn_unnormed = (time_term + radial_term + theta_term + phi_term) 

print ndotn_unnormed
B2 = sym.simplify(1/ndotn_unnormed)
B = sym.sqrt(B2)
print "The normalization constant is", B
print "In latex format:", sym.latex(sym.collect(B,a)) 

# Hint: To print expressions in latex format:
#latex(Integral(cos(x)**2, (x, 0, pi)))

# Hint: To enable pretty printing, write 
#init_printing(use_unicode=True)

# Hint: To change the value of a Symbol in an expression, use subs:
#expr.subs(x, 2)

# Hint: Good way to check if two expressions are equal:
#simplify(a - b)
# (In which case the result should be 0.)

# Hint: In order to avoid division problems, you can use the function rational:
#Rational(1, 2)
# Here you would get the rational 1/2.

# Here are some examples of expression substitution:
#>>> expr = x**y
#>>> expr
#x**y
#>>> expr = expr.subs(y, x**y)
#>>> expr
#x**(x**y)
#>>> expr = expr.subs(y, x**x)
#>>> expr
#x**(x**(x**x))

#>>> expr = x**3 + 4*x*y - z
#>>> expr.subs([(x, 2), (y, 4), (z, 0)])
#40

# Hint: Here is an example of how to 
#convert a SymPy expression to an expression that can be numerically evaluated:
#>>> import numpy 
#>>> a = numpy.arange(10) 
#>>> expr = sin(x)
#>>> f = lambdify(x, expr, "numpy") 
#>>> f(a) 
#[ 0.          0.84147098  0.90929743  0.14112001 -0.7568025  -0.95892427
# -0.2794155   0.6569866   0.98935825  0.41211849]

# Hint: collect() collects common powers of a term in an expression
#collect() is particularly useful in conjunction with the .coeff() method. 
#expr.coeff(x, n) gives the coefficient of x**n in expr