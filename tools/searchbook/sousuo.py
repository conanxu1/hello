import requests

from urllib import parse
from PyQt5.QtWidgets import  \
                        QApplication,\
                        QLabel,\
                        QWidget,\
                        QVBoxLayout,\
                        QLineEdit,\
                        QTextEdit,\
                        QFormLayout,\
                        QPushButton,\
                        QListView,\
                        QMessageBox\


from PyQt5.QtCore import Qt, QStringListModel  
import sys
from bs4 import BeautifulSoup
import re





daohangwangzhi=[
"http://shu.ziyuandi.cn/",
"https://zhuanlan.zhihu.com/p/479925339",
"https://www.shuxiangjia.cn",
"https://www.showbooks.net/",
"https://www.gpdf.net/",
"https://findbooks.eu.org/",
"https://yabook.blog/",
"https://www.iyd.wang/",
"https://www.jingdoo.net/all",
"https://yabook.org/",
"https://vikingcabin.com/",
"http://www.dushudaren.com/",
"https://orzbooks.com/",
"http://www.33file.com/",
"https://linbaoshuwu.com/",
"http://www.eybook.com/",
"https://book.519.best/",
"https://www.tianlangbooks.com/",
"https://ibooks.org.cn",
"https://dianshu.xyz",
"https://www.nmod.net/",
"https://sobooks.net/",
"https://www.soushu.vip/",
"https://www.youyidu.xyz/",
"https://www.yueduxiu.vip/",
"https://www.book123.info/",
"http://www.happydot.top",
"https://ebooklist.mobi/",
"https://xxurls.com/",
"http://www.78books.com/",
"https://aishuchao.com/",
"https://www.kdlbook.cn/",
"http://www.ziliaoh.com/",
"http://www.rejoiceblog.com/",
"https://http561856124.wordpress.com/",
"https://bookdown.org/",
"http://kgbook.com/",
"https://www.shuge.org/collections/",
"http://www.ixinqing.com/",
"https://pdfzj.com/",
"http://www.woaidu.cc/",
"https://zh.annas-archive.org/",
"",
"",
]


 
global shoulis
shoulis=[
"https://www.gpdf.net",
"https://www.duolaweb.com",
"https://www.ershu.org/",
"https://www.cnblogs.com/eson15/p/12653193.html",
"https://new.shuge.org/view/yue_jue_shu/",
"https://www.gokanla.com/6732.html",
"https://sobooks.net",
"https://www.jiumodiary.com",
"http://bks.thefuture.top/",
"http://gutenberg.org"
"https://www.shuzhou.cc/"
"https://www.iyd.wang/"
"https://book.zhishikoo.com/"
"https://downloadpdf.cn/index.html"
"http://www.pdfzj.cn/"
"http://shu.ziyuandi.cn/"
"https://www.jb51.net/books/"
"https://www.junengfan.cn/bd"
"https://www.chaonengso.com/"
"https://homeofpdf.com/"
"https://pdfshare.top"
]   
            


def f1(sou,shuming,wwiinn):


    res=requests.get(sou+"?q="+shuming) 
    wwiinn.te1.setPlainText(res.text)




# # # global urllis
# # # urllis=[]
# # # for ee in shoulis:
    # # # urllis.append(ee+"/index.php?m=vod-search")

 





# # # def souv1():
    # # # global yingpian 
    # # # global shou
    # # # global url
    # # # data = {"wd":yingpian,
            # # # "submit":"search"}
    # # # res = requests.post(url=url,data=data)

    
    # # # return res,shou,'utf-8'








# 可能要登陆

 
 
 



# # # # def wolongzy():
    # # # # global yingpian 
    # # # # url = "https://wolongzy.net/search.html"
    # # # # data = {"searchword":yingpian}
    # # # # res = requests.post(url=url,data=data)

    # # # # html=res.text
    # # # # return html,'utf-8'


# # # # # 请求网址:https://wolongzy.net/search.html?searchword=%E6%9F%AF%E5%8D%97



# # # # def bajieziyuan():
    # # # # global yingpian 
    # # # # url = "http://bajieziyuan.com/index.php?m=vod-search"
    # # # # data = {"wd":yingpian}
    # # # # res = requests.post(url=url,data=data)
    # # # # html=res.text
    # # # # return html,'utf-8'

# # # # def kuyuntv():
    # # # # global yingpian 
    # # # # url = "https://kuyun.tv/vod/search.html"
    # # # # data = {"wd":yingpian}
    # # # # res = requests.post(url=url,data=data)
    # # # # html=res.text
    # # # # return html,'utf-8'




 



# # # # def kuyunzy1():
    # # # # global yingpian 
    # # # # url = "http://www.kuyunzy1.com/search.asp"
    # # # # data = 'searchword='+parse.quote(yingpian.encode("gb2312")) 

    # # # # header = {'content-type':'application/x-www-form-urlencoded'}

 
    # # # # res = requests.post(url=url,data=data, headers=header)

    # # # # html=str(res.text.encode('latin1').decode('gb2312'))
    # # # # return html, "gb2312"

 
 
# 可能要浏览器登陆一下
























class Demo(QWidget):
    def __init__(self):
        global shoulis

        
        super().__init__()
        qfl = QFormLayout()
        self.setMinimumWidth(800)
        self.setMaximumWidth(5000)
        self.setMinimumHeight(600)
        self.setMaximumHeight(1000)
        
        
        
        
        self.le1 = QLineEdit()
        self.le1.setPlaceholderText("书名")
        qfl.addRow("Normal", self.le1)
        
        
        
        # self.bt1 = QPushButton('PyQt5 button', self)
        # self.bt1 .setToolTip('This is an example button')
        # self.bt1 .clicked.connect(self.on_click)
        # qfl.addRow("bt", self.bt1 )
        
        
        
        
        self.lsv1=QListView()
        #实例化列表模型，添加数据
        self.slm=QStringListModel()
        self.qList=shoulis
        self.slm.setStringList(self.qList)
        self.lsv1.setModel(self.slm)
        self.lsv1.clicked.connect(self.clicked)
        qfl.addRow("liebiao", self.lsv1 )
        
        self.te1=QTextEdit()
        qfl.addRow("反馈", self.te1)
        
        
        
        
        
        EDIT='''
        self.lsv2=QListView()
        #实例化列表模型，添加数据
        self.slm2=QStringListModel()
        self.qList2=["ok","ok"]
        self.slm2.setStringList(self.qList2)
        self.lsv2.setModel(self.slm2)
        self.lsv2.clicked.connect(self.clicked2)
        qfl.addRow("liebiao2", self.lsv2 )
        
        self.le2=QLineEdit()
        qfl.addRow("地址", self.le2)
        
        
        
        
        
        
        
        self.lsv3=QListView()
        #实例化列表模型，添加数据
        self.slm3=QStringListModel()
        self.qList3=["m3u8","mp4","refresh"]
        self.slm3.setStringList(self.qList3)
        self.lsv3.setModel(self.slm3)
        self.lsv3.clicked.connect(self.clicked3)
        qfl.addRow("提取", self.lsv3 )
        '''

        
        
        
        self.setLayout(qfl)
        
         
        
    def clicked(self,qModelIndex):
        #提示信息弹窗，你选择的信息
        idd=qModelIndex.row()
        global shou
        global url
        global yingpian
        global shoulis
        global urllis
        global piandan1
        global piandan2
        
        piandan1=[]
        piandan2=[]
        
        sou=shoulis[idd]
        # # # # url=urllis[idd]
        shuming=self.le1.text()
        
        
        try:
            
            turl="https://www.shuxiangjia.cn/search.php"
            if(sou==turl):
                print("get")
                f1(sou,shuming,self)
                
                
                

                
        
        except:
            print("err")
        
        
        
        aaa='''
        try:
            res,shou,ccode= souv1()
            soup = BeautifulSoup(res.content, "html.parser")
            lis1=soup.find_all(class_="xing_vb4")
            for ee in lis1:
                # print(ee.text)
                # print(shou+ee.a.get("href"))
                piandan1.append(ee.text)
                piandan2.append(shou+ee.a.get("href"))
                
                
            self.qList2=piandan1
            self.slm2.setStringList(self.qList2)
            # self.lsv2.setModel(self.slm2)
                
                
        except:
            print("\n++er++\n")
        
        
        
        
        print("\n+++++++++++++++++++++++++++++\n")
        # print(res.text)
        # print("\n+++++++++++++++++++++++++++++\n")
        # QMessageBox.information(self,'ListWidget','你选择了：'+self.qList[qModelIndex.row()])
        '''
        
        
    
    
  
    # # # # # # # # # def clicked2(self,qModelIndex):
        # # # # # # # # # global piandan2
        
  
        # # # # # # # # # idd2=qModelIndex.row()
        # # # # # # # # # res = requests.get(piandan2[ idd2])
        # # # # # # # # # self.te1.setPlainText(res.text)
        # # # # # # # # # self.le2.setText(piandan2[ idd2]) 
        
    
    # # # # # # # # # def clicked3(self,qModelIndex):
        # # # # # # # # # strr=self.te1.toPlainText()
        # # # # # # # # # shuchu=""
        # # # # # # # # # if self.qList3[qModelIndex.row()]=="mp4":
            # # # # # # # # # jieguo=re.findall('''(http.*?mp4)''', strr)
            # # # # # # # # # jieguosort = list(set(jieguo))
            # # # # # # # # # jieguosort.sort(key=jieguo.index)
            
            # # # # # # # # # for ee in jieguosort:
                # # # # # # # # # shuchu=shuchu+ee+"\n"
                
                
                
            # # # # # # # # # self.te1.setPlainText(shuchu)
        # # # # # # # # # elif self.qList3[qModelIndex.row()]=="m3u8":
            # # # # # # # # # jieguo=re.findall('''(http.*?m3u8)''', strr)
            # # # # # # # # # jieguosort = list(set(jieguo))
            # # # # # # # # # jieguosort.sort(key=jieguo.index)
            
            # # # # # # # # # for ee in jieguosort:
                # # # # # # # # # shuchu=shuchu+ee+"\n"
            # # # # # # # # # self.te1.setPlainText(shuchu)               
        # # # # # # # # # elif self.qList3[qModelIndex.row()]=="refresh":
            # # # # # # # # # try:
                # # # # # # # # # res = requests.get( self.le2.text()  )
                # # # # # # # # # self.te1.setPlainText(res.text)
            # # # # # # # # # except:
                # # # # # # # # # print("er")
           

 
        
        
        
        
 

if __name__=="__main__":
    app=QApplication(sys.argv)
    win=Demo()
    win.show()
    sys.exit(app.exec_())


 

# while(1):
    # yingpian=input("mingzi\n")
    # html,shou,ccode= zuidazy2()

# f=open('tt.html','w',encoding=ccode)
# f.write(html)
# f.close()


# for i in L:
# ...     if not i in T:
# ...         T.append(i)






 