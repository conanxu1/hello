from bs4 import BeautifulSoup
from bs4 import UnicodeDammit
import re
import gzip
import re
import http.cookiejar
import urllib.request
import urllib.parse
import socket
import os
import time





def downloadPage(url):
    h = urllib.request.urlopen(url)
    return h.read()

def downloadImg(content):
    pattern = '''(http.*?")'''
    m = re.compile(pattern)
    soup = BeautifulSoup(content, 'html.parser') 
    
    urls = soup.find_all("img")
    # print(soup)
    f=open("my.txt","a+",encoding="utf8")
    
    for e in urls:
    
        f.write(e["src"].replace("?","\n\n")+"\n\n")
    f.close()
        
    
    
    # for i, url in enumerate(urls):
        # print(url)
        
        
        # # urllib.urlretrieve(url, "%s.jpg" % (i, ))



for i in range(1,20):
    content = downloadPage("https://naixi.lofter.com/"+"?page="+str(i) )
    downloadImg(content)


# 


