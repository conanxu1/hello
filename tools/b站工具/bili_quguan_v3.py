from DrissionPage import ChromiumPage
from DrissionPage import SessionPage
from DrissionPage import WebPage
import json
from bs4 import BeautifulSoup
import time






def mget_data(uid,csrf):
    data={"fid" :  uid,
      "act" :   "2",
      "csrf":	csrf,
      }
    return data


 
f=open("name_id.txt","r",encoding="utf8")
w=f.readlines()
f.close()



urla="https://space.bilibili.com/485156027/fans/follow"
#用于登录



page = WebPage()
print(page.mode)

time.sleep(1)
page.get(urla)
input("登录......\n")
page.change_mode()

csrf=page.cookies(as_dict=True)["bili_jct"]
for ee in w:
    aaa=(ee.replace("\n","").split(";;;"))
    uid=str(aaa[0])
    data=mget_data(uid,csrf)
    page.post("https://api.bilibili.com/x/relation/modify",data=data)
    time.sleep(0.2)
    print(uid)


