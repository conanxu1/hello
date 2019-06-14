
from cvxopt import matrix, solvers
import numpy as np


Q = 2*matrix([[5, -.5], [-.5, 1]])
p = matrix([-2, -6.0],(2,1))




#dengshi
A = matrix([[4.0,5.0]],(1,2))
b = matrix([7.0])

#budengshi    liexiangliang
G = matrix([[1.0, -1],[8,5]])
h = matrix([2,2.0],(2,1))

print(A)

sol=solvers.qp(Q, p, G, h, A, b)
print(sol['x'])
print(sol['primal objective'])
print(sol)