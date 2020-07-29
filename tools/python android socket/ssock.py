# python 3.6
# 多层代理
# 功能：转发至socks5server

from socketserver import ThreadingMixIn, TCPServer, StreamRequestHandler
import socketserver # socketserver- 网络服务器的框架
import socket # socket- 低级网络接口
import struct # struct - 将字节解释为压缩二进制数据
import time # time - 时间访问和转换
import select # select - 等待I / O完成
import selectors # selectors - 高级I / O复用
import logging

logging.basicConfig(level=logging.DEBUG)



class ThreadingTCPServer(ThreadingMixIn, TCPServer):
    pass

class Encoder(StreamRequestHandler):

    def exchange_loop(self, client, remote):
        while True:

            # wait until client or remote is available for read
            r, w, e = select.select([client, remote], [], [])

            if client in r:
                data = client.recv(4096) # 请求数据包
                if remote.send(data) <= 0:
                    break

            if remote in r:
                data = remote.recv(4096) # 响应数据包
                if client.send(data) <= 0:
                    break


    def handle(self):
        logging.info('Accepting connection from %s:%s' % self.client_address)

        #: 来自客户端的连接
        local = self.connection
        #: 连接至SOCKS server
        remote = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        remote.connect(('127.0.0.1', 51234))
        self.exchange_loop(local, remote)



def main():
    server = ThreadingTCPServer(('127.0.0.1', 51234), Encoder)
    server.serve_forever()

if __name__ == '__main__':
    main() 
