import dryscrape
from bs4 import BeautifulSoup

dryscrape.start_xvfb()
session = dryscrape.Session()
my_url = 'https://space.bilibili.com/4856007/video'
session.visit(my_url) 
response =session.body() 
soup = BeautifulSoup(response,features="lxml")
soup.find(id="intro-text")
print(soup)


































