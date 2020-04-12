# pip install PyQtWebEngine
# setUrl setHtml
# QLabel, ,QFormLayout 
#winshell win32com


#专门为rsshub设置导入功能

import time
from PyQt5.QtCore import QT_VERSION_STR
from PyQt5.Qt import PYQT_VERSION_STR
from sip import SIP_VERSION_STR
import sys
from PyQt5.QtWidgets import QApplication,\
                            QMainWindow,\
                            QWidget,\
                            QHBoxLayout,\
                            QVBoxLayout,\
                            QLineEdit,\
                            QTextEdit,\
                            QPushButton,\
                            QListView,\
                            QMessageBox

from PyQt5 import QtWebEngineWidgets
                           
                         
from PyQt5.QtCore import QUrl,\
                         QStringListModel  
                         
from PyQt5.QtWebEngineWidgets import QWebEngineView
import os
import feedparser
import json
import socket
socket.setdefaulttimeout(25)

import ssl

 

# 2. 表示忽略未经核实的SSL证书认证

context = ssl._create_unverified_context()
 



global rssbiao
global dingyuebiao
global wenzhangbiao
global d

rssbiao=[]
dingyuebiao=[]
wenzhangbiao=[]

def fresh():
    global rssbiao
    global dingyuebiao
    rssbiao=[]
    dingyuebiao=[]
    f=open("rssbiao.txt","r",encoding="utf-8")
    w=f.readlines()
    f.close()
    
    tw=[]
    for ee in w:
        tw.append("".join(ee.split()))

    nw = [ee for ee in tw if ee != '']
    for ee in nw:
        tem=ee.split("==")
        print(tem)
        dingyuebiao.append(tem[0].replace(" ", ""))
        rssbiao.append(tem[1].replace(" ", ""))
    

class WebEngineView(QtWebEngineWidgets.QWebEngineView):
    def createWindow(self,QWebEnginePage_WebWindowType):
        page = WebEngineView(self)
        page.urlChanged.connect(self.on_url_changed)
        return page
    def on_url_changed(self,url):
        self.setUrl(url)


 
class MainWindow(QMainWindow):
    def __init__(self):
        global rssbiao

        super(MainWindow, self).__init__()
        self.setWindowTitle('加载外部网页的例子')
        self.setGeometry(5,30,1355,730)
        self.browser=WebEngineView()
        #加载外部的web界面
        self.browser.load(QUrl('https://blog.csdn.net/jia666666'))
        
        
        
        
        
        
        self.btn1=QPushButton('打开链接')
        self.btn1.clicked.connect(self.on_click_fresh)
        
        self.btn2=QPushButton('打开ssr')
        self.btn2.clicked.connect(self.on_click_ssr)
        
        self.btn3=QPushButton("更新列表")
        self.btn3.clicked.connect(self.on_click_freshdingyue)
        
        
        #订阅信息
        self.lsv1=QListView()
        self.slm=QStringListModel()
        self.qList=dingyuebiao
        self.slm.setStringList(self.qList)
        self.lsv1.setModel(self.slm)
        self.lsv1.clicked.connect(self.clicked)
        # 文章列表  
        self.lsv2=QListView()
        self.slm2=QStringListModel()
        self.qList2=["订阅号"]
        self.slm2.setStringList(self.qList2)
        self.lsv2.setModel(self.slm2)
        self.lsv2.clicked.connect(self.clicked2)
      
        self.le2=QLineEdit()
        
        self.te1=QTextEdit()
        self.te1.setPlainText("相关信息")
      
      
        
        
        widget = QWidget()
        self.setCentralWidget(widget)   # 建立的widget在窗体的中间位置
 
        self.hBox = QHBoxLayout()
        self.vBox= QVBoxLayout()
        
        self.vBox.addWidget(self.te1)
        self.vBox.addWidget(self.btn3)
        self.vBox.addWidget(self.lsv1)
        self.vBox.setStretchFactor(self.te1,1)
        self.vBox.setStretchFactor(self.lsv1,3)
        
        
        self.vBox2= QVBoxLayout()
        self.vBox2.addWidget(self.le2)
        self.vBox2.addWidget(self.btn1)
        self.vBox2.addWidget(self.btn2)
        self.vBox2.addWidget(self.lsv2)
  
  
              
        self.hBox.addLayout(self.vBox)
        self.hBox.addLayout(self.vBox2)
        self.hBox.addWidget(self.browser)
        self.hBox.setStretchFactor(self.vBox,1)
        self.hBox.setStretchFactor(self.vBox2,2)
        self.hBox.setStretchFactor(self.browser,7)
        
        
        widget.setLayout(self.hBox)
        
        
    def clicked(self,qModelIndex):
        
        global wenzhangbiao
        global d
        
        wenzhangbiao=[]
        
        #提示信息弹窗，你选择的信息
        idd=qModelIndex.row()
        
        d=feedparser.parse( rssbiao[idd]  )
        
        try:
            self.le2.setText( d.feed.link)
            strr=d.feed.title+"\n" + d.feed.subtitle+"\n"  
            self.te1.setPlainText(strr)
        except:
            self.te1.setPlainText("er!!")
        
        for ee in d.entries:
            wenzhangbiao.append( ee.title )
        
        self.qList2=wenzhangbiao
        self.slm2.setStringList(self.qList2)
        
        

    def clicked2(self,qModelIndex):
        global d
        
        idd=qModelIndex.row()
        
        tem=d.entries
        ttiid=0
        
        try:
            for ee in d.entries:
                if ttiid==idd:
                    tem=ee.links
                    uurl=tem[0].href
                ttiid=ttiid+1
           
            self.le2.setText(uurl)
            self.browser.setUrl(QUrl(uurl))
        except:
            self.le2.setText("error")
            ttid=0
            for ee in d.entries:
                if ttiid==idd:
                    tem=ee.links
                    print(tem[0])
                    if tem[0].href !="":
                        self.le2.setText(tem[0].href)
                    
                ttiid=ttiid+1
        
        
    def on_click_fresh(self):
        uurl=self.le2.text()
        try:
            self.browser.setUrl(QUrl(uurl))
            
        except:
            yuan=self.te1.toPlainText()
            print(yuan)
            self.te1.setPlainText(yuan+"\n无法打开\n")
      
      
    
    def on_click_ssr(self):
        res=os.popen("D:\小工具\ShadowsocksR-4.6.1\ShadowsocksR-dotnet4.0.exe")
    
    def on_click_freshdingyue(self):
        global dingyuebiao
        fresh()
        self.qList=dingyuebiao
        self.slm.setStringList(self.qList)

 
if __name__ == '__main__':
    fresh()
    app=QApplication(sys.argv)
    win=MainWindow()
    win.show()
    print("Qt5 Version Number is: {0}".format(QT_VERSION_STR))
    print("PyQt5 Version is: {}".format(PYQT_VERSION_STR))
    print("Sip Version is: {}".format(SIP_VERSION_STR))


    app.exit(app.exec_())
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
# class MainWindow(QMainWindow):
  # def __init__(self, *args, **kwargs):
    # super().__init__(*args, **kwargs)
    # self.setWindowTitle('My Browser')
    # self.showMaximized() 
    # self.webview = WebEngineView()
    # self.webview.load(QUrl("https://www.bing.com"))
    # self.setCentralWidget(self.webview)
 
# class WebEngineView(QWebEngineView):
  # windowList = []
  # def createWindow(self, QWebEnginePage_WebWindowType):
    # new_webview =  WebEngineView()
    # new_window = MainWindow()
    # new_window.setCentralWidget(new_webview)
    # self.windowList.append(new_window)
    # return new_webview
 
# if __name__ == "__main__":
  # app = QApplication(sys.argv)
  # window = MainWindow()
  # window.show()
  # sys.exit(app.exec_()) 
  
  
  
  
# import sys
 
 
# from PyQt5.QtCore import QUrl,QApplication
# from PyQt5.QtWebKit import QWebView
# from PyQt5.QtGui import QGridLayout, QLineEdit, QWidget
 
# class UrlInput(QLineEdit):
    # def __init__(self, browser):
        # super(UrlInput, self).__init__()
        # self.browser = browser
        # # add event listener on "enter" pressed
        # self.returnPressed.connect(self._return_pressed)
 
    # def _return_pressed(self):
        # url = QUrl(self.text())
        # # load url into browser frame
        # browser.load(url)
 
# if __name__ == "__main__":
    # app = QApplication(sys.argv)
 
    # # create grid layout
    # grid = QGridLayout()
    # browser = QWebView()
    # url_input = UrlInput(browser)
    # # url_input at row 1 column 0 of our grid
    # grid.addWidget(url_input, 1, 0)
    # # browser frame at row 2 column 0 of our grid
    # grid.addWidget(browser, 2, 0)
 
    # # main app window
    # main_frame = QWidget()
    # main_frame.setLayout(grid)
    # main_frame.show()
 
    # # close app when user closes window
    # sys.exit(app.exec_())
