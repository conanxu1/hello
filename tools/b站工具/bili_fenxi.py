import time
# from pyvirtualdisplay import Display
import json

from DrissionPage import WebPage
from DrissionPage import ChromiumOptions
from DrissionPage import ChromiumPage
from DrissionPage.easy_set import set_paths
from DrissionPage.easy_set import set_headless

import os



#co=ChromiumOptions().set_paths(browser_path="/usr/bin/chromium")
#page = ChromiumPage(addr_driver_opts=co)
#page = ChromiumPage(addr_driver_opts=co)
#os.system("pkill chromium")



# # # # # display = Display(visible=0, size=(1280, 768))
# # # # # display.start()



port=9333
os.system("nohup chromium --no-sandbox  --remote-allow-origins=\"*\" --remote-debugging-port="+str(port)+"  >  log.txt  2>&1   &  ")

##co=ChromiumOptions().set_paths(local_port=port)
##page = WebPage(addr_driver_opts=co)



set_paths(local_port=port)
page = WebPage()




 




f=open("name_id_read.txt","r",encoding="utf-8")
w=f.readlines()
f.close()












page.get("http://bilibili.com")

input("====\n")
    

# uul="https://space.bilibili.com/"+str(aaa[0])+"/video"


# page.get(uul)


for ee in w:
    aaa=ee.split(";;;")
    print(aaa[1])
    print(aaa[0])
    
    
    
    uul="https://space.bilibili.com/"+str(aaa[0])+"/video"
    page.get(uul)
    
    
    
    time.sleep(5)
    
    uul="https://api.bilibili.com/x/space/wbi/arc/search?"\
        +"mid="+str(aaa[0])\
        +"&ps="+"10"\
        +"&pn=1&order=pubdate"

    
    page.get(uul)
    # # print(page.html)
    
    # js=json.loads(page.html.text)
    
    js=page.json
    
    print(js["data"]["list"]["vlist"])
    
    # # # # # # # print(page.html)
    
    
    # # # # # # # pp=page.ele('@class=content').eles("tag:ul")
    
    
    
    # # # # # # ##选择第一个隐藏的！！
    
    
    # # # # # # pp1=page.ele('@class=content')
    # # # # # # print(page.html)
    
    # # # # # # pp=page.ele('@class=content').eles("tag:ul")
    # # # # # # print(pp[1].html)
    
    # # # # # # # # pp=pp.ele("@class=list")
    
    
    # # # # # # for uuii in pp[1].eles("tag:li"):
        # # # # # # print(uuii)
        # # # # # # print("\n")
        

    time.sleep(10)
    

