#coding:utf-8

import socket, sys  
      
# Step1: 输入host和port信息  
host = "192.168.1.103"  #input('please input host name: ')  
textport ="51500"           #input('please input textport: ')  
      
# Step2: 创建socket对象  
s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)  
try:  
    port = int(textport)  
except ValueError:  
    port = socket.getservbyname(textport, 'udp')  
       
# Step3: 打开socket连接     
s.connect((host, port))  
      
# Step4: 发送数据  
print("Enter data to transmit: ")  
data = sys.stdin.readline().strip()  
s.sendall(data.encode("utf8"))  
      
# Step5: 接收服务器发过来的数据  
print("Looking for replies; press Ctrl-C or Ctrl-Break to stop") 
while 1:  
    buf = s.recv(2048).decode()  
    if not len(buf):  
        break  
    sys.stdout.write(buf)  
