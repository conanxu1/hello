import re

f=open("new.txt","r",encoding="utf-8")
w=f.read()
f.close()


res=re.findall('''(www.bilibili.com/video.*?")''',w)

res=set(res)

f=open("out.txt","a+",encoding="utf-8")


for ee in res:
    f.write(ee.replace('''"''',"\n").replace('''www.''',"annie http://www."))
 


f.close()
