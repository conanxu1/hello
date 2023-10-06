
import socket
import socks

from bs4 import BeautifulSoup

import urllib3    # 导入urllib3模块
headers={'User-Agent' : "Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/49.0.2623.221 Safari/537.36 SE 2.X MetaSr 1.0 "}
    
url = "https://www.socks-proxy.net/"
url = "https://spys.one/en/socks-proxy-list/"
url3 = "https://freeproxyupdate.com/socks5-proxy"
 


http = urllib3.PoolManager()   # 创建连接池管理对象


socks.set_default_proxy(socks.SOCKS5, "127.0.0.1",10808)
socket.socket = socks.socksocket

r = http.request('GET', url,headers=headers)    # 发送GET请求

print(r.status,r.data.decode('utf-8')) 

soup=BeautifulSoup(r.data.decode('utf-8'),'lxml')


lst = soup.find_all('tr')

f=open("ipip.txt","w",encoding="utf-8")
for ee in lst:
    lstgg=ee.find_all('td')
    print("==========")
    
    
    
    count=0
    for gg in lstgg:
    
        print(gg)
        if count < 2:
            f.write(gg.text.strip())
            f.write(":")
            count += 1   
    f.write("\n")
    
    
    # for gg in lstgg:
        # print(gg.text)

f.close()

















