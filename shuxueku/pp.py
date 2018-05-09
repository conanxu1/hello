import matplotlib.pyplot as plt
import numpy as np

x = np.arange(9)
y = np.sin(x)
z = np.cos(x)
# marker数据点样式，linewidth线宽，linestyle线型样式，color颜色
plt.plot(x, z)
plt.title("matplotlib")
plt.xlabel("height")
plt.ylabel("width")
# 设置图例
plt.legend(["Y","Z"], loc="upper right")
plt.grid(True)



登录
hupozhao
关注
树莓派最新raspbian系统换国内源
转载 2017年10月10日 22:33:07 阅读 4809
树莓派新版系统更换了专门优化过的桌面环境PIXEL，正好手头有个闲置的TF卡决定刷上新版系统玩玩。下载刷系统过程很多教程页很简单。插卡，上电开机，释放卡上的剩余空间都很正常，因为树莓派官方源访问很慢下一步就是换成国内源，以前一直在用中科大的源，于是开始按照以前的方法修改/etc/apt/sources.list换源。可是换完发现还会从官方源更新，直到找到了/etc/apt/sources.d/raspi.list。

网上教程一般只更换/etc/apt/sources.list里面的源地址，当按照之前的方法更换/etc/apt/sources.d/raspi.list里面的网址之后，发现更新失败了...

直到从网上查找资料之后发现正确的方法：

 

1.修改之前，先备份下配置文件。

 

cp /etc/apt/sources.list
cp /etc/apt/sources.d/raspi.list
 


2.修改/etc/apt/sources.list或者直接修改原文件，把原有的配置全部注释掉（使用#注释）。

#deb http://mirrordirector.raspbian.org/raspbian/ jessie main
