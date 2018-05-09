from sympy.diffgeom import *
import sympy as sp
import numpy as np
from sympy.simplify import simplify
from sympy.diffgeom.rn import R3

#设定初始的符号 参数 坐标系抽象函数

V=sp.Symbol('V')
d0=sp.Symbol('d0')
d1=sp.Symbol('d1')
d2=sp.Symbol('d2')
u=sp.Symbol('u')
b=sp.Symbol('b')
q=sp.Function('q')


n = 3
M = Manifold('M', n)
P = Patch('P', M)

coord = CoordSystem('coord', P, ['x%s'%i for i in range(n)])
x = coord.coord_functions()
dx = coord.base_oneforms()
ex = coord.base_vectors()


fm=np.mat([[d0*(-x[0]-x[1]*x[2]+V)],[-d1*q(x[1])],[d2*(x[0]*x[1]-b*x[2])]])

gm=np.mat([[0],[d1],[0]])




f=np.mat(ex)*fm
g=np.mat(ex)*gm


f=f[0,0]
g=g[0,0]
h=x[2]


print(f)
print(g)
print(h)

jie=input('jie')
jie=int(jie)

t=h
for i in range(jie-1):
    t=LieDerivative(f,t).doit()

t2=LieDerivative(g,t).doit()

t1=LieDerivative(f,t).doit()




print(t1)
print(t2)

ans=t1+u*t2


sp.pprint(simplify(ans))


