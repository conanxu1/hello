
#coding:utf-8 
import requests
from lxml import etree
import urllib.request as request
import urllib.parse
from bs4 import BeautifulSoup
import time
import requests
import re
import binascii
import csv
import xlsxwriter
import difflib
 
import smtplib  
from email.mime.text import MIMEText  
from email.header import Header  
  
# 第三方 SMTP 服务  
mail_host="smtp.163.com"  #设置服务器  
mail_user="xuwenhan_2018@163.com"    #用户名  
mail_pass="MIPGOVZGMKLAWMLZ"   #口令,QQ邮箱是输入授权码，在qq邮箱设置 里用验证过的手机发送短信获得，不含空格  




sender = 'xuwenhan_2018@163.com'    # 发件人邮箱(最好写全, 不然会失败)  
receivers = [
					'872144755@qq.com',
					]  # 接收邮件，可设置为你的QQ邮箱或者其他邮箱  
 


def sendEmail(title,content):

    message = MIMEText(content, 'html', 'utf-8')  # 内容, 格式, 编码
    message['From'] = "{}".format(sender)
    message['To'] = ",".join(receivers)
    message['Subject'] = title

    try:
        smtpObj = smtplib.SMTP_SSL(mail_host, 465)  # 启用SSL发信, 端口一般是465
        smtpObj.login(mail_user, mail_pass)  # 登录验证
        smtpObj.sendmail(sender, receivers, message.as_string())  # 发送
        print("mail has been send successfully.")
    except smtplib.SMTPException as e:
        print(e)

import socket

def get_host_ip(): 
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(('8.8.8.8', 80))
        ip = s.getsockname()[0]
    finally:
        s.close()
    return ip

def pub_ip():
	try:
		res=requests.get("https://ifconfig.co/ip")
		print(res.text)
		
		print(res)
		return res.text
	
	except:
		print("er")
	
	finally:
		print("ok")
 

pub_ip()



sendEmail(get_host_ip(),"ok")
