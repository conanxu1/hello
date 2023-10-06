# import urllib.request

# proxy_host = '221.211.62.4'
# proxy_port = 1111
# proxy_handler = urllib.request.ProxyHandler({'socks5': f'{proxy_host}:{proxy_port}'})
# opener = urllib.request.build_opener(proxy_handler)
# response = opener.open('https://ifconfig.me/ip')
# print(response.read().decode())



# # import socket
# # import socks
# # import requests

# # # socks.set_default_proxy(socks.SOCKS5, "85.239.55.203",49572)
# # # socket.socket = socks.socksocket
# # print(requests.get('http://cip.cc').text) 





from tqdm import tqdm, trange

from urllib3.contrib.socks import SOCKSProxyManager
import socket
import socks
import urllib3    # 导入urllib3模块
headers={'User-Agent' : "Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/49.0.2623.221 Safari/537.36 SE 2.X MetaSr 1.0 "}
    
#uurl = "https://cip.cc"
uurl = "http://www.httpbin.org/ip"
# uurl = "https://ip.900cha.com/"




# # # def shi(url,port):
    # # # http = urllib3.PoolManager(timeout=urllib3.Timeout(connect=0.5, read=1.0))   # 创建连接池管理对象

    # # # socks.set_default_proxy(socks.SOCKS5 , url , port)
    # # # socket.socket = socks.socksocket
    
    # # # print(url)
    # # # r = http.request('GET', uurl,headers=headers)    # 发送GET请求
    # # # print(r.status,r.data.decode('utf-8')) 



def shi(url,port):
    proxy = SOCKSProxyManager("socks4://"+url+":"+str(port)+"/")
    proxy = SOCKSProxyManager("socks5://"+url+":"+str(port)+"/")
    r = proxy.request('GET', uurl ,headers=headers,timeout=1)
    
    
    print(r.status,r.data.decode('utf-8')) 




f=open("ipip.txt","r",encoding="utf-8")
ww=f.readlines()
f.close()



f=open("workip.txt","w",encoding="utf-8")
f.close()


for ee in tqdm(ww):
    kk=ee.split(":")
    f=open("workip.txt","a+",encoding="utf-8")

        
    try:
        print(kk[0],int(kk[1]))
        shi(kk[0],int(kk[1]))
        f.write(kk[0])
        f.write(":")
        f.write(kk[1])
        f.write("\n")
        
        
        
    except BaseException as e:
        a=0
    f.close()





# # # 
# # # 
# # # 







