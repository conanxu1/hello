import os
import numpy as np


f=open('log.txt','r',encoding='utf-8')

w=f.readlines()

g=0
n=np.shape(w);
print(n[0])



x=[]
for i in range(0,int(n/3)):

