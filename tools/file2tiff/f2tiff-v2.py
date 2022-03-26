import numpy as np
from PIL import Image
import os
import datetime
import math
from skimage import io
import time
import struct
import sys




#按字节写入
# f=open(file,"ab")
# f.write(struct.pack('B', num))
#f.close()
	
def mywrite(file,num):
	file.write(struct.pack('B', num))
    

 


#时间戳
def timestamp():
	return str(datetime.datetime.now()).replace(" ","-").replace(":","-").replace(".","-")







#文件转图片
#infosize:文件信息长度默认200
#file_mat:图片矩阵
#filesize:文件大小
#num1:所需矩阵元素个数，奇数时元素个数为num1+num2
#num2:是否为奇数
#tonji:记录所转换的字节个数

#ps:数字的utf8码为二位数，逗号的utf8码也为二位数，故可直接删除"0x"

def file2mat(urlname,png_wid,m_ts):
	global infosize
	
	file_mat=np.array(np.uint16(np.zeros((png_wid+1,png_wid)) ))

	filesize= os.path.getsize(urlname)
	
	#计算所需元素个数
	num1=math.floor(filesize/2)
	print(num1)
	num2=filesize%2
	print(num2)

    
    
    
	numsize=png_wid*png_wid
	 
	numk=math.ceil( filesize/(2*numsize)   )
	print("\nnumk:")
	print(numk)
	print("=================")
	f=open(urlname,"rb")
	
	tongji=0
	
	for k in range(1,numk+1):
	
		strr=str(k)+","+m_ts
		sbm=bytes(strr, encoding = "utf8")  
		jj=0
		for ee in sbm:
			file_mat[  png_wid  ,jj]=ee
			jj=jj+1
        
#sbm=bytes(strr, encoding = "utf8") 直接输出为b'文字'。若for ee in sbm:输出为数字
#hex(数字)="0x11"，int("0x11"，16)=数字

		print(k)
		
		
		for i in range(png_wid):
			for j in range(png_wid):
				
				# print(str(i)+","+str(j)+","+str(k)+"\n")
				w=f.read(2)

				
				c=0
				ll=0
				for x in w:
					#print(x)
					# c=c+int(hex(x),16)*(2**(ll) )
					c=c+x*(2**(ll) )
					ll=ll+8
				file_mat[i,j]=c
				#若num2=1此循环只执行一次
                
				

				
				
				tongji=tongji+1
				if tongji==num1+num2:
					break
			if tongji==num1+num2:
					break
				
		im=Image.fromarray(file_mat)
		im.save(m_ts+"="+str(k)+"="+".tiff")
		
		if tongji==num1+num2:
			break
	f.close()
	
	
	
	
	
	
	#文件信息 info.tiff------------------------------
	strr=str(num1)+","+str(num2)+","+str(numk)+","+str(png_wid)+","+urlname+","+m_ts
	sbm=bytes(strr, encoding = "utf8")  
 
	file_mat=np.array(np.uint16(np.zeros((1,infosize)) ))
	j=0
	for ee in sbm:

		file_mat[0,j]=ee
		j=j+1
		

	im=Image.fromarray(file_mat)
	im.save(m_ts+"=info"+".tiff")
	
	infoname=(m_ts+"=info"+".tiff")
	return infoname





#返回图片最后一行的信息
def endstr(img,png_wid):
	bb=''
	for jk in range(png_wid):
		addnew=hex(int(img[png_wid,jk]))
		bb=bb+addnew
		

	nstr=bb.rstrip("0x0").replace("0x","")

	try:
		byarray=bytearray.fromhex(nstr)
	except:
		byarray=bytearray.fromhex(nstr+"0")
	infostr=byarray.decode("utf8")             
	
	
	return infostr


	
	
	
#图转文件
def png2file(pngname):
	global infosize
    
    
    
    #读图--------------
	img = io.imread(pngname)
	# # # print(img)
	
    
	
    
    
    
    #解析文件信息--------
	bb=''
	for j in range(infosize):
		addnew=hex(int(img[0,j]))
		bb=bb+addnew

    #此处数字和逗号utf8码为两位数
	nstr=bb.rstrip("0x0").replace("0x","")
   
	try:
		byarray=bytearray.fromhex(nstr)
	except:
		byarray=bytearray.fromhex(nstr+"0")

    
    
	infostr=byarray.decode("utf8")
	infolist=infostr.split(",")
	
	num1		=int( infolist[0])
	num2		=int( infolist[1])
	numk		=int( infolist[2])
	png_wid	=int( infolist[3])
	filename	=infolist[4] 
	m_ts			=infolist[5]
	
	print(filename)
	print(numk)
	print("=================")
#############################################


#修改图片文件名
	dirs=os.listdir("./")
	for file in dirs:
		if file.endswith(".tiff"):
				img=io.imread(file)
				try:
					infosstr=endstr(img,png_wid)
					infolist=infosstr.split(",")
					if infolist[1]==	m_ts	:
						os.rename(file,m_ts+"="+infolist[0]+"="+".tiff")
				except:
					oooo=1



#恢复的文件名
	bkname="bk_"+filename
	f=open(bkname,"wb")
	
	
 

	# # # # k<numk range(1,1)为空
	count=0
	for k in range(1,numk):
		print(k)
		img = io.imread(m_ts+"="+str(k)+"="+".tiff")
		for i in range( png_wid):
			for j in range( png_wid):
 
				tmpbyte2= math.floor(  int(img[i,j])/2**8 )  
				tmpbyte1= math.floor(  int(img[i,j])%2**8 )  
 

				mywrite(f,tmpbyte1)
				mywrite(f,tmpbyte2)
	
	n2num1=num1-(numk-1)*png_wid*png_wid



	numi=math.floor(   n2num1/png_wid  )
    #多出的行
	numj=math.floor(   n2num1%png_wid  )
	#多出的列
 
 
#k=1时也执行
	for k in range(numk,numk+1):
		print(k)    
		img = io.imread(m_ts+"="+str(k)+"="+".tiff")
		for i in range( numi):
			for j in range( png_wid):
 
			
				tmpbyte2= math.floor(  int(img[i,j])/2**8 )  
				tmpbyte1= math.floor(  int(img[i,j])%2**8 )  
 
				mywrite(f,tmpbyte1)
				mywrite(f,tmpbyte2)
				
		for i in range( numi,numi+1):
			for j in range( numj):
				

				tmpbyte2= math.floor(  int(img[i,j])/2**8 )  
				tmpbyte1= math.floor(  int(img[i,j])%2**8 )  

				mywrite(f,tmpbyte1)
				mywrite(f,tmpbyte2)


		
		if num2>0:
			tmpbyte=math.floor(  int(img[numi,numj])%2**8 )  
			mywrite(f,tmpbyte)


 

	f.close()	
	print("ok")    
		
	
	
	
    
print("python new.py filename 3000 backup"+"\n"+"python new.py filename.tiff recover"+"最小值200")    
    
	
    
    
    
    
    
	
global infosize
infosize =200
	
	






try:
	urlname=sys.argv[1]
	png_wid=int(sys.argv[2])
	m_ts=timestamp()
except:

	ooo=1

try:
	if sys.argv[3]=="backup":
		infoname=file2mat(urlname,png_wid,m_ts)
		print("ok")
		
except  Exception as e:
	print("error:")
	print(e)

if sys.argv[2]=="recover":
	print("ok")
	png2file(sys.argv[1])
				
	
	
	
	
	
	
#================================================	
	
	
# def info2png(urlname):
	# filesize= os.path.getsize(urlname)
	
# i,j,k
def  next(i,j,k,nn):
	if (i+1)*(j+1)==nn*nn:
		i=0
		j=0
		k=k+1
	elif (j+1)==nn:
		i=i+1		
		j=0
	else:
		j=j+1
	return i,j,k
		
def next_ij(i,j,nn):
	if j==nn-1:
		i=i+1
		j=0
	else:
		j=j+1
	return i,j

# def myread(m_ts,i,j,k)
	
	
	
	
	
	
	
	# # for i in range(nn):
		# # for j in range(nn):
			# # w=f.read(2)
			# # print(w)
			# # if w ==b'':
				# # return file_mat
			# # c=0
			
			
			# # ll=8
			# # for x in w:
				# # c=c+int(hex(x),16)*(2**(ll) )
				# # ll=ll-8
			# # print(c)
			# # file_mat[i,j]=65535
			# # im.save(timestamp()+"="+str(k)+".tiff")

			
# # file_mat=file2mat(f)
# # im=Image.fromarray(file_mat)





# # # # # urlname="LU.txt"
# # # # # filesize= os.path.getsize(urlname)


# # # # # print( filesize )



# # # # # f=open(urlname,"rb")

# # # # # for itr in range(filesize   ):


# # # # # w=f.read(1)
# # # # # print(w)

# # # # # # num=int.from_bytes(w, byteorder='big')
# # # # # # print(num)
# # # # # for x in w:
	# # # # # print(hex(x))


# # # # # w=f.read(1)
# # # # # print(w)
# # # # # for x in w:
	# # # # # print(hex(x))
# # # # # w=f.read(1)
# # # # # print(w)
# # # # # for x in w:
	# # # # # print(int(hex(x),16))


# #im =Image.fromarray(a, mode='I;16')
# im.save("test.tiff")

# print(timestamp())





# # per_size=nn*nn*2

# # file_mat=np.array(np.uint16(np.ones((nn,nn)) ))

# f=open(urlname,"rb")
# filesize= os.path.getsize(urlname)
# f.close()

# from PIL import Image

# img = Image.new('L', (200 , 50), (0))
# img.save("test.png", "PNG", transparency=0)

