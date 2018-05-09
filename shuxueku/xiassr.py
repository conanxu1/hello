# -*- coding:utf-8 -*-
from selenium  import  webdriver

option = webdriver.ChromeOptions()
option.add_argument('--headless')
driver = webdriver.Chrome(chrome_options=option)
# driver = webdriver.Chrome()
# driver = webdriver.PhantomJS()
driver.get('https://www.baidu.com/')
print('打开浏览器')
print(driver.title)
driver.find_element_by_id('kw').send_keys('测试')
print('关闭')
driver.quit()
print('测试完成')
