

from bs4 import BeautifulSoup

import urllib3    # 导入urllib3模块
headers={'User-Agent' : "Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/49.0.2623.221 Safari/537.36 SE 2.X MetaSr 1.0 "}
    

#url2 = "https://proxy-tools.com/proxy/socks5?page="



   # 创建连接池管理对象
# # socks.set_default_proxy(socks.SOCKS5, "85.239.55.203",49572)
# # socket.socket = socks.socksocket
  



def zhuijia(url):

    http = urllib3.PoolManager(timeout=8.0)
    
    
    try:
        r = http.request('GET', url,headers=headers)
        soup=BeautifulSoup(r.data.decode('utf-8'),'lxml') 
        lst = soup.find_all('tr')

        # 发送GET请求
    except BaseException as e:
        lst=[]
    
    # print(r.status,r.data.decode('utf-8')) 

    


    print(len(lst))

    f=open("ipip.txt","a+",encoding="utf-8")
    for ee in lst:
        lstgg=ee.find_all('td')
        
        
        
        count=0
        for gg in lstgg:
        
            
            if count < 2:
                f.write(gg.text.strip())
                f.write(":")
                count += 1   
        f.write("\n")
        
        
        # for gg in lstgg:
            # print(gg.text)

    f.close()


f=open("ipip.txt","w",encoding="utf-8")
f.close()



def zhuijia2(url,n):
    for i in range(1,n):
        zhuijia(url+str(i))
        print("=+"+str(i))



url1 = "https://www.proxy-list.download/SOCKS5"


url2="https://www.freeproxy.world/?type=socks5&page="
url3="https://freeproxyupdate.com/socks5-proxy"



print("1")
zhuijia(url1)
print("2")
zhuijia(url3)
print("3")
zhuijia2(url2,3)












