'''直接使用bs4分析不保存页面，并研究如何分析微盘里面的文件夹

正则匹配最大页数
<a href="?page=1"
设计递归函数搜索目录
le:表征是第几级目录
全局变量技术。控制le   变量
'''


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
 

def last(url):

    url1=url+'1'
    data=urllib.request.urlopen(url1).read()
    soup=BeautifulSoup(data,'html.parser')
    result = re.findall(re.compile(r'href="\?page=(.*?)"'), str(soup))
    return int(result[len(result)-2])





def folder(urlt,le):
    
#得到对应地址的html    
    data=urllib.request.urlopen(urlt).read()
    soup=BeautifulSoup(data,'html.parser')
#打开记录文件
    time.sleep(3)
    
    reObj=re.compile(r'sort_name_pic([\s\S]*?)class="(.*?)"([\s\S]*?)class="sort_name_intro"([\s\S]*?)href="(.*?)"([\s\S]*?)title="(.*?)"')
    result = re.findall(reObj, str(soup))
    f=open(file,'a+',encoding='utf-8')
    f.write('page='+str(i)+'\n')
    f.close()


    
    for  a in result:
        
        print(i)
        
        kind=a[1]
        urll=a[4]
        name=a[6]
        f=open(file,'a+',encoding='utf-8')


        
        if kind!='vd_icon32_v2 vd_folder':

        
            f.write(name)
            f.write('\n')
            f.close()
            

            
        else:
           
            f.write('lev='+str(le))
            
 
            le=le+1
            f.write('+++++')

            f.write(name)
            
            f.write('\n')
            f.close()

            folder(urll,le)

            f=open(file,'a+',encoding='utf-8')
            f.write(']]]]]]]lev='+str(le)+'\n\n')
            f.close()


        le=0   
    return 1
    










def page(url,i):
      
    le=0
    
    folder(url+str(i),le)
    print('+++++++++++++++++++')


     







'''初始化参数'''

url='http://vdisk.weibo.com/u/3237600867?page='
file='好书推荐-阴阳俞.txt'



'''求最大页数'''
#n=last(url)
#print(n)
n=200




'''递归调用输出目录'''


for i in range(3,n+1):
    #print(i)
    page(url,i)
    print('++++++++++++++++++++++++++++++++++++++++++++')











#f=open(file,'w+',encoding='utf-8') 
#f.close()
'''
   
    f1.write(str(i)+'.................\n\n')
    f1.close()  
    print(i)
    time.sleep(5)
    return 1

'''

#    soup = BeautifulSoup(writer,"lxml")
#   for url in soup.findAll(href=re.compile("^qb:(.*)"))

   

'''
reObj= re.compile(r'vd_icon32_v2 vd_folder([\s\S]*?)title="(.*?)"')
    result=re.findall(reObj,str(soup))
    if len(result)>0:
        for a in result:
            print(a[1])
    
'''    

