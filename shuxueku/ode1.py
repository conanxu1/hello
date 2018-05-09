import numpy as np 
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from numpy import linspace 


miu=0.1
b=0.03
c=1

def sat(g):
    if np.abs(g)<=1:
        return g
    else:
        return np.sign(g)
def r(t):
    return np.sin(t/3)
def u(e1,e2):
    return -(2*np.abs(e2)+3)*sat((e1+e2)/miu)


def f(w, t): 
    # 给出位置矢量w，和三个参数p, r, b计算出
    # dx/dt, dy/dt, dz/dt的值


    x1,x2,e1,e2 = w
    # 直接与lorenz的计算公式对应 
    dx1=x2
    dx2=-np.sin(x1)-b*x2+c*u(e1,e2)
    #de1=x2-np.cos(t/3)/3
    de1=e2
    de2=-np.sin(x1)-b*x2+c*u(e1,e2)+np.sin(t/3)/9
    
    return np.array([dx1,dx2,de1,de2]) 

times = linspace(0.0001, 10, 10**5) 



sol = odeint(f, (1.5,0,1.5,-1/3), times) 


t=times
x = np.sin(times/3) 


plt.plot(t,x,t,sol[:,0]) 
plt.show()


