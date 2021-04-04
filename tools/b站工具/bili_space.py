from selenium import webdriver
from selenium.webdriver.common.action_chains import ActionChains
from selenium.webdriver.common.keys import Keys
from bs4 import BeautifulSoup

import time



def bili_space():
	




	
	
	action.key_down(Keys.SPACE).key_up(Keys.SPACE).perform()
	action.key_down(Keys.SPACE).key_up(Keys.SPACE).perform()
	action.key_down(Keys.SPACE).key_up(Keys.SPACE).perform()
	page = driver.page_source
	soup=BeautifulSoup(page,'html.parser')	
	array=soup.findAll('a',attrs={"class":"title"})
	
	for ee in array:
		jieguo_id.append(ee.get("href"))
		jieguo_name.append(ee.text)
		

	down_btn =  driver.find_element_by_xpath("//li[@class='be-pager-next']")
	down_btn.click()

	time.sleep(3)








profile = webdriver.FirefoxProfile()
profile.set_preference("general.useragent.override", "Mozilla/5.0 (Windows NT 6.1; rv:22.0) Gecko/20130405 Firefox/22.0")
fireFoxOptions = webdriver.FirefoxOptions()
# fireFoxOptions.headless = True
driver = webdriver.Firefox(profile,options=fireFoxOptions)
action = ActionChains(driver)


# urla=input("我的关注\n")


jieguo_id=[]
jieguo_name=[]
urla="https://space.bilibili.com/485156027/fans/follow"
driver.get(urla)
input("......")	

for i in range(106):
	print(i)
	try:
		bili_space()
	except:
		pass

f=open("name_id.txt","a",encoding="utf-8")
print(jieguo_id)
print(len(jieguo_id))
for i in range(len(jieguo_id)):
	f.write(jieguo_id[i]+";;;"+jieguo_name[i]+"\n")
f.close()


driver.close()
driver.quit()
	





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # =
# # # # chromeOptions = webdriver.ChromeOptions()
# # # # # chromeOptions.add_argument("--headless")

# # # # chromeOptions.add_argument('--no-sandbox')
# # # # chromeOptions.add_argument('--disable-dev-shm-usage')


# # # # driver = webdriver.Chrome(options=chromeOptions)
# # # # action = ActionChains(driver)


# # # # bili_space()


# # # # driver.close()
# # # # driver.quit()