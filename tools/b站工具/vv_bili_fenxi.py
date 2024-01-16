from DrissionPage import WebPage
import time
from pyvirtualdisplay import Display
from DrissionPage.easy_set import set_paths


from DrissionPage import ChromiumOptions
from DrissionPage.easy_set import set_headless 
set_headless(True)

from DrissionPage import ChromiumPage, ChromiumOptions



display = Display(visible=0, size=(1280, 768))
display.start()

do1 = ChromiumOptions().set_paths(local_port=9111)

# 创建多个页面对象
page1 = ChromiumPage(addr_driver_opts=do1)





 
co = ChromiumOptions()
co.set_argument('--incognito')
co.set_argument('--no-sandbox') 
 
set_paths(browser_path="/usr/bin/chromium") #一般linux安装的google浏览器默认都在这个目录




f=open("name_id_read.txt","r",encoding="utf-8")
w=f.readlines()
f.close()


for ee in w:
    aaa=ee.split(";;;")
    print(aaa[1])
    

