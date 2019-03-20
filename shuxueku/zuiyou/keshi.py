import os
import numpy as np
import matplotlib.pyplot as plt



from pyvirtualdisplay import Display

#是否虚拟输出
if 1:
	display = Display(visible=0, size=(800, 600))
	display.start()








f=open('log.txt','r',encoding='utf-8')


w=f.readlines()
f.close()
g=0
n=np.shape(w);
print(n[0])

ttemp=w[0].split(',')

print(ttemp)
print(np.shape(ttemp)[0])
dim=np.shape(ttemp)[0]-1


x=[]
t=[]

for i in range(0,int(n[0])):
	temp=w[i].split(',')
	x.append(float(temp[0].rstrip(",")))
	
	
	
	t.append(float(temp[dim].lstrip("??").rstrip("\n")))
	print(t[i])
	print("\n\n")
	
	
	

plt.plot(t,x)
plt.show()

plt.savefig("3.jpg")
plt.close()
