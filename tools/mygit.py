import re
from selenium import webdriver


browser = webdriver.Chrome()

def fenxi(browser,uurl):
    
    browser.get(uurl)
    js = "return document.documentElement.outerHTML"
    html = browser.execute_script(js)

    res=re.findall('''(www.bilibili.com/video.*?")''',html)
    res=set(res)

    f=open("out.txt","a+",encoding="utf-8")
    for ee in res:
    f.close()




turl="https://space.bilibili.com/396521874/video"


for i in range(2,5):
    
    fenxi(browser,turl+"?page="+str(i))



