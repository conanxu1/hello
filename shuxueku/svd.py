#coding:utf-8  
'''
求矩阵奇异值分解小应用



array对象和mat对象数积矢量积
shape求行列'''




 
from numpy import *  
  
def loadData():  
    return [[39.926,	41.68],
			[52.688,65.6860]
				
			]  
  
data=loadData()  
  
u,sigma,vt=linalg.svd(data)  
  
sig=zeros([u.shape[0],vt.shape[0]]) 


print(sigma);
for i in range(1,sigma.shape[0]+1):
	sig[i-1,i-1]=sigma[i-1]
	
	
print(mat(u)*mat(sig)*mat(vt))