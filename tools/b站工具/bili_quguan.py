
import os
import time
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as ec
from bs4 import BeautifulSoup

from selenium.webdriver.common.action_chains import ActionChains
from selenium.webdriver.common.keys import Keys

import requests
import re

import json
import random
from urllib import request
from urllib import parse
from urllib import request
import sys
import urllib


def page_bili(uurl):
	time.sleep(2)
	state=0
	
	
	
	driver.get(uurl)
	html = driver.page_source
	soup=BeautifulSoup(html,'html.parser')
	array=soup.findAll('img')
	f=open("bili_pic.sh","a",encoding="utf-8")
	for ee in	array:
		try:	
			f.write("aria2c -c -x 16 -s 16 http:"+ee.get("src")+"\n")
		except:
			pass
		try:	
			f.write("aria2c -c -x 16 -s 16 http:"+ee.get("data-src")+"\n")
		except:
			pass
			
			
			
	f.close()


	# # # # # # res=re.findall("(//.*?jpg)",html)
	# # # # # # res2=re.findall("(//.*?png)",html)



	
	# # # # # # for ee in res:
		
	# # # # # # for ee in res2:
		# # # # # # f.write("aria2c -c -x 16 -s 16 http:"+ee+"\n")

 





#
#def bili():
#	
#
#	jieguo=[]
#
#	#################
#
#	 
#    
#	print("开始批量取关")
#	
#	pp=100
#	
#	driver.get("http://baidu.com")
#	action = ActionChains(driver)
#		
#	time.sleep(300)
#		
#    
#    
#    
#    
#    
#    
#    
#	jieguo=[]
#	
#	
##############################################	
#	
#	
#	oo=1
#	for i in range(pp):
#		print(pp)
#		try:
#		
#			page = driver.page_source
#			soup=BeautifulSoup(page,'html.parser')
#			
#			
#			if bili_opt==1:
#				array=soup.findAll('h2',attrs={"class":"article-title"})
#			if bili_opt==2:
#				array=soup.findAll('a',attrs={"class":"title"})
#			
#			# print(array)
#			
#			for ee in array:
#			
#				try:
#					
#					
#					
#					if bili_opt==1:
#						jilu=ee.find("a").get("href")
#						
#						jieguo.append("http:"+jilu+'\n')
#					if bili_opt==2:
#						jilu=ee.get("href")
#						print(jilu)
#					
#						jieguo.append("annie -p http:"+jilu+'\n')
#			
#				except:
#					print(ee)
#					print("\n")
#
#				
#			print(i)	
#			print("\n")
#				
#			action.key_down(Keys.SPACE).key_up(Keys.SPACE).perform()
#			time.sleep(5)
#			action.key_down(Keys.SPACE).key_up(Keys.SPACE).perform()	
#			time.sleep(5)
#			if oo<pp:
#				down_btn =  driver.find_element_by_xpath("//li[@class='be-dropdown-item']")
#		
#				time.sleep(5)
#				down_btn.click()
#				oo+=1
#		except:
#			print("....")
#			
#			
#
#
#	if bili_opt==1:
#	
#	
#	
#
#		for ee in jieguo:
#			print(ee)
#			page_bili(ee)
#
#		
#		
#		
#		
#		
#	
#	if bili_opt==2:
#		f=open("bili_vid.sh","w",encoding="utf-8")
#		f.write(" ".join(jieguo))
#		f.close()
#	
##os.chdir("D:\mygit\hello\tools\b站工具")
#print(os.getcwd())	
	
profile = webdriver.FirefoxProfile()
profile.set_preference("general.useragent.override", "Mozilla/5.0 (Windows NT 6.1; rv:22.0) Gecko/20130405 Firefox/22.0")
fireFoxOptions = webdriver.FirefoxOptions()
#fireFoxOptions.headless = True
#driver = webdriver.Firefox(profile,options=fireFoxOptions)
driver = webdriver.Firefox()
	


driver.get("https://space.bilibili.com/485156027/fans/follow")
action = ActionChains(driver)

input("开始取关")

for i in range(200):
    print(i)
    driver.get("https://space.bilibili.com/485156027/fans/follow")
    time.sleep(0.2)
    driver.execute_script('''$(".be-dropdown-item:contains('取消关注')").click()''')
    time.sleep(0.2)










driver.close()
driver.quit()
	











# from requests_html import HTMLSession

# session = HTMLSession()


# uurl=input("dizhi")
# r = session.get(uurl)
 
# print(r.html.html)