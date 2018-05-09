import tushare as ts
import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.api as sm
import matplotlib.pyplot as  plt

df=ts.get_k_data(code='600010',start='2017-01-01',end='2018-01-01')


dta=df['open']
dta=np.array(dta,dtype=np.float)


print(dta)
arma_70=sm.tsa.ARMA(dta,(7,0)).fit()


resid=arma_70.resid
print(resid)


plt.clf()
plt.plot(resid)
plt.plot(dta)


plt.show()


