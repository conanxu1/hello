from DrissionPage import ChromiumPage
from DrissionPage import SessionPage
from DrissionPage import WebPage


import json


from bs4 import BeautifulSoup

import time




'''

order_type=attention&
gaia_source=main_web&
web_location=333.999&
w_rid=40bcfa8cc539ec7cb444cf52aab3502a&
wts=1702777802
'''

jieguo_id=[]
jieguo_name=[]
 


def bili_space_n(page,i):
    uu="https://api.bilibili.com/x/relation/followings?"\
        +"vmid=485156027"\
        +"&ps=50&order=desc"\
        +"&pn="+str(i)
    
   
    page.get(uu)
     
    
    
    js=json.loads(page.html)
    
    
    f=open("name_id.txt","a+",encoding="utf-8")

    for ee in (js["data"]["list"]):
        # print(ee)

       
        print(ee["mid"])
        print(ee["uname"])
        jieguo_id.append(ee["mid"])
        jieguo_name.append(ee["uname"])
        f.write(str(ee["mid"])+";;;"+ee["uname"]+"\n")
    f.close() 
        



 





urla="https://space.bilibili.com/485156027/fans/follow"




page = WebPage()
print(page.mode)

time.sleep(3)
page.get(urla)
input("......")	
page.change_mode()
print(page.mode)
print(page.title)
time.sleep(2)






f=open("name_id.txt","w",encoding="utf-8")
f.close()



for i in range(70,71):
    
    print("++start_num\n")
    print(i)
    
    
    try:
        bili_space_n(page,i)
        print("++end_num\n")
        print(i)

        time.sleep(5)
    except:
        pass


 
# page.quit()




