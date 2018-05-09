#coding:utf-8  
'''
求矩阵奇异值分解小应用



array对象和mat对象数积矢量积
shape求行列'''




 
from numpy import *  
  
def loadData():  
    return [[0,-1.6,0.6],
			[0,1.2,0.8],
			[0,0,0],
			[0,0,0]
				
			]  
  
data=loadData()  
  
u,sigma,vt=linalg.svd(data)  
  
sig=zeros([u.shape[0],vt.shape[0]]) 



for i in range(1,sigma.shape[0]+1):
	sig[i-1,i-1]=sigma[i-1]
	
	
print(mat(u)*mat(sig)*mat(vt))