import os
import time
import tushare as ts
import pandas as pd

from youjian import sendEmail

def check(code, low, high):
	df = ts.get_realtime_quotes(code)
	e = df[['code','name','price','time']]
	p = df[u'price']
	
	print(e) 
	
	
	if float(p[0]) > low and float(p[0]) < high:
		return True
	else :
		return False
	
liebiao=[
		['sh', 2800, 3100],
		['600486',39,60]
		]	
	
while True:
	
	
	
	
	for ee in liebiao:
		if not check(ee[0], ee[1], ee[2]):
			print(1)
			sendEmail("gupiao",ee[0]+"\n")
		
		
		# # exit()
	time.sleep(1500)
