import tushare as ts  
import pandas as pd  
import numpy as np  
import scipy.stats as scs #科学计算  
import matplotlib.pyplot as plt #绘图  
import statsmodels.api as sm #统计运算 
  
a=ts.get_hist_data('600848',start='2016-11-01',end='2017-01-02')  
a=a['close']  
a.name = '600848'  
  
b=ts.get_hist_data('000002',start='2016-11-01',end='2017-01-02')  
b=b['close']  
b.name='000002'  
  
c=ts.get_hist_data('002285',start='2016-11-01',end='2017-01-02')  
c=c['close']  
c.name='002285'  

print(c)
  
  
d=pd.DataFrame([a,b,c])  
##转置  
data=d.T  
##回报率  
returns = np.log(data / data.shift(1))  
##年化收益率  
returns.mean()*252  
##计算协方差矩阵  
returns.cov()*252  
  
##计算股票个数  
noa=len(data.T)  


print(noa)


##随机生成初始化权重  
weights = np.random.random(noa)  
##计算百分比  
weights /= np.sum(weights)  
weights  
  
  
##下面通过一次蒙特卡洛模拟，产生大量随机的权重向量，并记录随机组合的预期收益和方差。  
port_returns = []  
  
port_variance = []  
  
for p in range(4000):  
    weights = np.random.random(noa)  
    weights /=np.sum(weights)  
    port_returns.append(np.sum(returns.mean()*252*weights))  
    port_variance.append(np.sqrt(np.dot(weights.T, np.dot(returns.cov()*252, weights))))  
##因为要开更号，所以乘两次weight  
##dot就是点乘   
port_returns = np.array(port_returns)  
port_variance = np.array(port_variance)  
  
#无风险利率设定为4%  
risk_free = 0.04  
plt.figure(figsize = (8,4))  
plt.scatter(port_variance, port_returns, c=(port_returns-risk_free)/port_variance, marker = 'o')  
plt.grid(True)  
plt.xlabel('excepted volatility')  
plt.ylabel('expected return')  
plt.colorbar(label = 'Sharpe ratio')  
plt.show()  
##投资组合优化1——sharpe最大  
  
def statistics(weights):  
    weights = np.array(weights)  
    port_returns = np.sum(returns.mean()*weights)*252  
    port_variance = np.sqrt(np.dot(weights.T, np.dot(returns.cov()*252,weights)))  
    return np.array([port_returns, port_variance, port_returns/port_variance])  
  
#最优化投资组合的推导是一个约束最优化问题  
import scipy.optimize as sco  
  
  
  
  
#最小化夏普指数的负值  
def min_sharpe(weights):  
    return -statistics(weights)[2]  
  
#约束是所有参数(权重)的总和为1。这可以用minimize函数的约定表达如下  
cons = ({'type':'eq', 'fun':lambda x: np.sum(x)-1})  
  
#我们还将参数值(权重)限制在0和1之间。这些值以多个元组组成的一个元组形式提供给最小化函数  
bnds = tuple((0,1) for x in range(noa))  
  
  
opts = sco.minimize(min_sharpe, noa*[1./noa,], method = 'SLSQP', bounds = bnds, constraints = cons)  
  
opts  
  
  
##sharpe最大的组合3个统计数据分别为：  
  
  
#预期收益率、预期波动率、最优夏普指数  
  
statistics(opts['x']).round(3)  
  
  
##通过方差最小来选出最优投资组合。  
  
  
def min_variance(weights):  
    return statistics(weights)[1]  
optv = sco.minimize(min_variance, noa*[1./noa,],method = 'SLSQP', bounds = bnds, constraints = cons)  
  
optv  
  
##方差最小的最优组合权重向量及组合的统计数据分别为：  
optv['x'].round(3)  
  
#得到的预期收益率、波动率和夏普指数  
statistics(optv['x']).round(3)
