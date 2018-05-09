from sympy.diffgeom import Manifold, Patch, CoordSystem
from sympy.diffgeom import TensorProduct as TP
from sympy.diffgeom import metric_to_Riemann_components as Riemann
from sympy.diffgeom import metric_to_Christoffel_2nd as Christoffel
from sympy import Symbol, Function, latex, sin

n = 2
M = Manifold('M', n)
P = Patch('P', M)

coord = CoordSystem('coord', P, ['x%s'%i for i in range(n)])
x = coord.coord_functions()
dx = coord.base_oneforms()
g = [[1, 0], [0, sin(x[0])**2]]
metric = sum([g[i][j]*TP(dx[i], dx[j]) for i in range(n) for j in range(n)])

C = Christoffel(metric)
R = Riemann(metric)