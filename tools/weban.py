from bs4 import BeautifulSoup
from bs4 import UnicodeDammit
import re
import gzip
import http.cookiejar
import urllib.request
import urllib.parse
import socket
import os
import time

def page(url,i,tmp):

    # urllib.request.urlretrieve(url+str(i),tmp)
    
    
    conn = urllib.request.urlopen(url+str(i)) 
    res=conn.read()

    # f=open(tmp,'r',encoding='utf-8')

    # res=f.read()
        # # reObj= re.compile('<span class="pic">(\r*)(\s*)<img src="(.*)"')
        # # result=reObj.findall(res) 
    # f.close()
    
    
    soup=BeautifulSoup(res,"html.parser")
 
    print(soup.get_text())
 
 
 
 
 
        # flag=0
        # while(flag==0):
            # try:
                # print('.')
                # urllib.request.urlretrieve(url+str(i),tmp)
                # flag=1
            # except:
                # time.sleep(20)
                # print('wait\n')
         
        # f=open(tmp,'r',encoding='utf-8')  
        # time.sleep(8)

        # res=f.read()
        # reObj= re.compile('<span class="pic">(\r*)(\s*)<img src="(.*)"')
        # result=reObj.findall(res) 
        # f.close()
        # time.sleep(8)





 

        # s=1
        
        # for a in result:
                # try:
                        # conn = urllib.request.urlopen(a[2])  
                # except:
                        # print('error')
                # f1 = open(str(i)+'-'+str(s)+'.jpg','wb')  
                # f1.write(conn.read())  
                # f1.close()  
                # print(str(i)+'.'+str(s)+'\n')
                # time.sleep(10)
                # s=s+1
        

        

        # print(i)
        # time.sleep(20)
        # f=open(tmp,'w+')
        # f.close()
        # return 1





tmp='tmp.html'


url='http://www.gyhj.org/thread/49/'

page(url,1,tmp)












#    soup = BeautifulSoup(writer,"lxml")
#   for url in soup.findAll(href=re.compile("^qb:(.*)"))






