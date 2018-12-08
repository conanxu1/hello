import os
import numpy as np


#from pyvirtualdisplay import Display
# display = Display(visible=0, size=(800, 600))
# display.start()


import matplotlib.pyplot as plt







f=open('log.txt','r',encoding='utf-8')


w=f.readlines()
f.close()
g=0
n=np.shape(w);
print(n[0])



x=[]
t=[]

for i in range(0,int(n[0]/3)):
	x.append(float(w[i*3].split(',')[0]))
	t.append(float(w[i*3].split(',')[1]))


plt.plot(t,x)
plt.show()

plt.savefig("3.jpg")
