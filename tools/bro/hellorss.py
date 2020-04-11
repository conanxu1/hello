import requests
from urllib import parse
from PyQt5.QtWidgets import  QApplication,QLabel,QWidget,QVBoxLayout,QLineEdit,QTextEdit,QFormLayout,QPushButton,QListView,QMessageBox
from PyQt5.QtCore import Qt, QStringListModel  
import sys
from bs4 import BeautifulSoup
import re
from graphviz import Digraph
import json
import socket
socket.setdefaulttimeout(25)



# 获取所有节点中最多子节点的叶节点
def getMaxLeafs(myTree):
    numLeaf = len(myTree.keys())
    for key, value in myTree.items():
        if isinstance(value, dict):
            sum_numLeaf = getMaxLeafs(value)
            if sum_numLeaf > numLeaf:
                numLeaf = sum_numLeaf
    return numLeaf


def plot_model(tree, name):
    g = Digraph("G", filename=name, format='png', strict=False)
    first_label = list(tree.keys())[0]
    g.node("0", first_label)
    _sub_plot(g, tree, "0")
    leafs = str(getMaxLeafs(tree) // 10)
    g.attr(rankdir='LR', ranksep=leafs)
    g.view()


root = "0"


def _sub_plot(g, tree, inc):
    global root

    first_label = list(tree.keys())[0]
    ts = tree[first_label]
    for i in ts.keys():
        if isinstance(tree[first_label][i], dict):
            root = str(int(root) + 1)
            g.node(root, list(tree[first_label][i].keys())[0])
            g.edge(inc, root, str(i))
            _sub_plot(g, tree[first_label][i], root)
        else:
            root = str(int(root) + 1)
            g.node(root, tree[first_label][i])
            g.edge(inc, root, str(i))

 

 



import feedparser
d=feedparser.parse("https://63c9b619.ngrok.io/bilibili/bangumi/media/9192")

print(d.feed )

print(d.feed.title )
print(d.feed.subtitle )
print(d.feed.link )


# for ee in d.entries : 
    # print(ee.title)
    # print(ee.summary_detail.value)
    # print(ee.id)
    
    # print("\n")


print(d.entries)

  
# for ee in d.entries: 
    # print("+++++++++++\n")
    
    # print(ee)
    # print(json.dumps(ee, indent=7))
    # for ff in ee: 
        # print(ff)
        # print("\n")
    # print("\n")
    # print("\n")
    # print("\n")
             






# # qita='''
# # https://kuyun.tv
# # "http://www.1886zy.net",
# # http://niuniuzy.com/index.php?m=vod-search

# # http://ivi.bupt.edu.cn/
# # https://www.mb700.com/vod/search.html

# # https://www.solezy.com/index.php/vod/search.html

# # '''



# # global shoulis
# # shoulis=[
# # "http://zuidazy2.com",
# # "http://zy.ataoju.com",
# # "http://www.765zy.com",
# # "http://www.156zy.me",
# # "http://www.666zy.com",
# # "https://www.subo988.com",
# # "http://www.mahuazy.net",
# # "http://www.haozy.cc",
# # "http://www.123ku.com",
# # "http://www.88zyw.net",
# # "http://265zy.cc",
# # "http://www.zuixinzy.net",
# # "http://kubozy.co",
# # "http://zy.itono.cn",
# # "http://kankanzy.com",
# # ]   
            


# # global urllis


# # urllis=[]
# # for ee in shoulis:
    # # urllis.append(ee+"/index.php?m=vod-search")

 





# # def souv1():
    # # global yingpian 
    # # global shou
    # # global url
    # # data = {"wd":yingpian,
            # # "submit":"search"}
    # # res = requests.post(url=url,data=data)

    
    # # return res,shou,'utf-8'








# # # 可能要登陆

 
 
 



# # def wolongzy():
    # # global yingpian 
    # # url = "https://wolongzy.net/search.html"
    # # data = {"searchword":yingpian}
    # # res = requests.post(url=url,data=data)

    # # html=res.text
    # # return html,'utf-8'


# # # 请求网址:https://wolongzy.net/search.html?searchword=%E6%9F%AF%E5%8D%97



# # def bajieziyuan():
    # # global yingpian 
    # # url = "http://bajieziyuan.com/index.php?m=vod-search"
    # # data = {"wd":yingpian}
    # # res = requests.post(url=url,data=data)
    # # html=res.text
    # # return html,'utf-8'

# # def kuyuntv():
    # # global yingpian 
    # # url = "https://kuyun.tv/vod/search.html"
    # # data = {"wd":yingpian}
    # # res = requests.post(url=url,data=data)
    # # html=res.text
    # # return html,'utf-8'




 



# # def kuyunzy1():
    # # global yingpian 
    # # url = "http://www.kuyunzy1.com/search.asp"
    # # data = 'searchword='+parse.quote(yingpian.encode("gb2312")) 

    # # header = {'content-type':'application/x-www-form-urlencoded'}

 
    # # res = requests.post(url=url,data=data, headers=header)

    # # html=str(res.text.encode('latin1').decode('gb2312'))
    # # return html, "gb2312"

 
 
# # # 可能要浏览器登陆一下
























# # class Demo(QWidget):
    # # def __init__(self):
        # # global shoulis

        
        # # super().__init__()
        # # qfl = QFormLayout()
        
        # # self.le1 = QLineEdit()
        # # self.le1.setPlaceholderText("影片")
        # # qfl.addRow("Normal", self.le1)
        
        
        
        # # # self.bt1 = QPushButton('PyQt5 button', self)
        # # # self.bt1 .setToolTip('This is an example button')
        # # # self.bt1 .clicked.connect(self.on_click)
        # # # qfl.addRow("bt", self.bt1 )
        
        
        
        
        # # self.lsv1=QListView()

        # # #实例化列表模型，添加数据
        # # self.slm=QStringListModel()
        # # self.qList=shoulis
        # # self.slm.setStringList(self.qList)
        # # self.lsv1.setModel(self.slm)
        # # self.lsv1.clicked.connect(self.clicked)
        # # qfl.addRow("liebiao", self.lsv1 )
        
        
        
        
        # # self.lsv2=QListView()

        # # #实例化列表模型，添加数据
        # # self.slm2=QStringListModel()
        # # self.qList2=["ok","ok"]
        # # self.slm2.setStringList(self.qList2)
        # # self.lsv2.setModel(self.slm2)
        # # self.lsv2.clicked.connect(self.clicked2)
        # # qfl.addRow("liebiao2", self.lsv2 )
        
        # # self.le2=QLineEdit()
        # # qfl.addRow("地址", self.le2)
        
        
        # # self.te1=QTextEdit()
        # # qfl.addRow("页面", self.te1)
        
        
        
        
        # # self.lsv3=QListView()

        # # #实例化列表模型，添加数据
        # # self.slm3=QStringListModel()
        # # self.qList3=["m3u8","mp4","refresh"]
        # # self.slm3.setStringList(self.qList3)
        # # self.lsv3.setModel(self.slm3)
        # # self.lsv3.clicked.connect(self.clicked3)
        # # qfl.addRow("提取", self.lsv3 )
    

        
        
        
        # # self.setLayout(qfl)
        
         
        
    # # def clicked(self,qModelIndex):
        # # #提示信息弹窗，你选择的信息
        # # idd=qModelIndex.row()
        # # global shou
        # # global url
        # # global yingpian
        # # global shoulis
        # # global urllis
        # # global piandan1
        # # global piandan2
        
        # # piandan1=[]
        # # piandan2=[]
        
        # # shou=shoulis[idd]
        # # url=urllis[idd]
        # # yingpian=self.le1.text()
        # # try:
            # # res,shou,ccode= souv1()
            # # soup = BeautifulSoup(res.content, "html.parser")
            # # lis1=soup.find_all(class_="xing_vb4")
            # # for ee in lis1:
                # # # print(ee.text)
                # # # print(shou+ee.a.get("href"))
                # # piandan1.append(ee.text)
                # # piandan2.append(shou+ee.a.get("href"))
                
                
            # # self.qList2=piandan1
            # # self.slm2.setStringList(self.qList2)
            # # # self.lsv2.setModel(self.slm2)
                
                
        # # except:
            # # print("\n++er++\n")
        
        # # print("\n+++++++++++++++++++++++++++++\n")
        # # # print(res.text)
        # # # print("\n+++++++++++++++++++++++++++++\n")
        
        
        
        # # # QMessageBox.information(self,'ListWidget','你选择了：'+self.qList[qModelIndex.row()])
 
    
    # # def clicked2(self,qModelIndex):
        # # global piandan2
        
  
        # # idd2=qModelIndex.row()
        # # res = requests.get(piandan2[ idd2])
        # # self.te1.setPlainText(res.text)
        # # self.le2.setText(piandan2[ idd2]) 
        
    
    # # def clicked3(self,qModelIndex):
        # # strr=self.te1.toPlainText()
        # # shuchu=""
        # # if self.qList3[qModelIndex.row()]=="mp4":
            # # jieguo=re.findall('''(http.*?mp4)''', strr)

            # # for ee in jieguo:
                # # shuchu=shuchu+ee+"\n"
            # # self.te1.setPlainText(shuchu)
        # # elif self.qList3[qModelIndex.row()]=="m3u8":
            # # jieguo=re.findall('''(http.*?m3u8)''', strr)
            # # for ee in jieguo:
                # # shuchu=shuchu+ee+"\n"
            # # self.te1.setPlainText(shuchu)               
        # # elif self.qList3[qModelIndex.row()]=="refresh":
            # # try:
                # # res = requests.get( self.le2.text()  )
                # # self.te1.setPlainText(res.text)
            # # except:
                # # print("er")
           


        
        
        
        
 

# # if __name__=="__main__":
    # # app=QApplication(sys.argv)
    # # win=Demo()
    # # win.show()
    # # sys.exit(app.exec_())


 

# # # while(1):
    # # # yingpian=input("mingzi\n")
    # # # html,shou,ccode= zuidazy2()

# # # f=open('tt.html','w',encoding=ccode)
# # # f.write(html)
# # # f.close()









 