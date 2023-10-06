
#tsocks proxychains  privoxy

url="https://proxylist.geonode.com/api/proxy-list?limit=500&sort_by=lastChecked&sort_type=desc&speed=fast&protocols=socks5"
    
    
import datetime  
from bs4 import BeautifulSoup
import urllib3    # 导入urllib3模块
headers={"User-Agent" : "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:109.0) Gecko/20100101 Firefox/118.0",
    "Referer":"https://geonode.com/",
    

}






#url2 = "https://proxy-tools.com/proxy/socks5?page="



   # 创建连接池管理对象
# # socks.set_default_proxy(socks.SOCKS5, "85.239.55.203",49572)
# # socket.socket = socks.socksocket


http = urllib3.PoolManager(timeout=8.0)
    
r = http.request('GET', url,headers=headers)



soup=BeautifulSoup(r.data.decode('utf-8'),'lxml')

import json
json1 = json.loads(soup.text)
print(json1.keys())


# jsso = sorted(json1["data"], key=itemgetter('lastChecked'))


f=open("ipip.txt","w",encoding="utf-8")
f.close()




for ee in (json1["data"]):
    
    
    f=open("ipip.txt","a+",encoding="utf-8")
    
    print("\n\n")
    print(ee.keys())
    
    print(ee["ip"])
    print(ee["port"])
    f.write(ee["ip"]+":"+ee["port"]+":"+"\n")
    f.close()
    
    print(ee["isp"])
    print(ee["city"])
    print(ee["latency"])
    
    
    
    timestamp=(ee["lastChecked"])
    dt_object = datetime.datetime.fromtimestamp(timestamp)
    formatted_time = dt_object.strftime('%Y-%m-%d %H:%M:%S')
    print(formatted_time)
    
    


