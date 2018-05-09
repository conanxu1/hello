from sympy import *

x=Symbol('x')
y=Symbol('y')
z=Symbol('z')

print(derive_by_array([x*sin(y),x*cos(y)],[x,y]))

