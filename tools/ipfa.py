
import socket
import fcntl
import struct
 
import os
from youjian import sendEmail
process = os.popen('sudo ifconfig') # return file
output = process.read()
process.close()
sendEmail("ip地址",output)


#
#def get_ip_address(ifname):
#    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
#    print(s)
#    print(s.fileno())
#    s.connect(('baidu.com', 80))
#    ip = s.getsockname()[0]
#    print(s.fileno())
#
#
#
#
#    #socket.inet_ntoa(
#    fcntl.ioctl(s.fileno(), 0x8915, struct.pack('256s', bytes(ifname[:15].encode("utf-8")))[20:24])
#     # SIOCGIFADDR
#    return 0
#
#
#
#get_ip_address('wlan0')
#
#
#
