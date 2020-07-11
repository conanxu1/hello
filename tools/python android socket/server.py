
#coding:utf-8

import socket, traceback  
      
host = '192.168.0.115' # Bind to all interfaces   
port = 1234  


#ifconfig
      
# Step1: 创建socket对象  
s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)  
      
# Step2: 设置socket选项(可选)  
s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)  
      
# Step3: 绑定到某一个端口  
s.bind((host, port))  
      
# Step4: 监听该端口上的连接  
while 1:  
    try:  
        message, address = s.recvfrom(8192)
        
        print("Got data from ", address,"1")  
        print(message.decode("utf8"))
        ans=input('Your reply')
        print(ans)


        s.sendto(ans.encode("utf8") ,address)



        
    except(KeyboardInterrupt, SystemExit):  
        print("raise")  
        raise
    except :  
        print("traceback")  
        traceback.print_exc()  
