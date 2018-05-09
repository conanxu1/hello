import tushare as ts
import matplotlib.pyplot as plt
import  matplotlib as mpl
import mpl_finance as mpf
import numpy as np
import pandas as pd

lra=ts.get_rrr()


print(lra[:1])


c=lra.iloc[:,2]
print(c.values)
n=np.shape(c.values)
print(n[0])
t=np.linspace(1,n[0],n[0])
plt.plot(t,c.values)
plt.show()
