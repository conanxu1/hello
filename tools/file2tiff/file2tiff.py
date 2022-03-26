import numpy as np
from PIL import Image
import os
import datetime
import math
from skimage import io
import time
import struct
import sys

def mywrite(file,num):

	# f=open(file,"ab")
	# f.write(struct.pack('B', num))
	# f.close()
	
	file.write(struct.pack('B', num))

	
	
	
	
	

	# newhex=hexnum.replace("0x","")
	# print(newhex)
	# byarray=bytearray.fromhex(str(newhex))
	# print(byarray)
	# for ee in byarray:
		# print(ee)
		# file.write(ee)



def timestamp():
	return str(datetime.datetime.now()).replace(" ","-").replace(":","-").replace(".","-")

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



def file2mat(urlname,png_wid,m_ts):
	global infosize
	
	file_mat=np.array(np.uint16(np.zeros((png_wid+1,png_wid)) ))
	i=0
	j=0
	
	filesize= os.path.getsize(urlname)
	
	
	num1=math.floor(filesize/2)
	print(num1)
	num2=filesize%2
	print(num2)
	
	numsize=png_wid*png_wid
	 
	numk=math.ceil( filesize/(2*numsize)   )
	print(numk)
	print("=================")
	f=open(urlname,"rb")
	
	tongji=0
	
	for k in range(1,numk+1):
	
		strr=str(k)+","+m_ts
		sbm=bytes(strr, encoding = "utf8")  
		jj=0
		for ee in sbm:
			# # # # print(ee)
			file_mat[  png_wid  ,jj]=int(hex(ee),16)
			jj=jj+1
		
	
	
		print(k)
		i=0
		j=0
		
		for i in range(png_wid):
			for j in range(png_wid):
				
				# print(str(i)+","+str(j)+","+str(k)+"\n")
				w=f.read(2)

				
				c=0
				ll=0
				for x in w:
					c=c+int(hex(x),16)*(2**(ll) )
					ll=ll+8
				file_mat[i,j]=c
				
				

				
				
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
	
	
	
	
	
	
	
	strr=str(num1)+","+str(num2)+","+str(numk)+","+str(png_wid)+","+urlname+","+m_ts
	sbm=bytes(strr, encoding = "utf8")  
	# # # print(sbm)
	file_mat=np.array(np.uint16(np.zeros((1,infosize)) ))
	j=0
	for ee in sbm:
		# # # # print(ee)
		file_mat[0,j]=int(hex(ee),16)
		j=j+1
		
	# # # print(file_mat)
	im=Image.fromarray(file_mat)
	im.save(m_ts+"=info"+".tiff")
	
	infoname=(m_ts+"=info"+".tiff")
	return infoname



def endstr(img,png_wid):
	bb=''
	for jk in range(png_wid):
		addnew=hex(int(img[png_wid,jk]))
		bb=bb+addnew
		
	nstr=bb.replace("0x","").rstrip("0")
	try:
		byarray=bytearray.fromhex(nstr)
	except:
		byarray=bytearray.fromhex(nstr+"0")
	infostr=byarray.decode("utf8")             
	
	
	return infostr


	
	
	

def png2file(pngname):
	global infosize
	img = io.imread(pngname)
	# # # print(img)
	
	bb=''
	for j in range(infosize):
 
		addnew=hex(int(img[0,j]))
		bb=bb+addnew

	print(bb)

	
    #nstr=bb.replace("0x","").rstrip("0")
	
	nstr=bb.replace("0x","").rstrip("0")
	print(nstr)
	
    
	try:
		byarray=bytearray.fromhex(nstr)
	except:
		byarray=bytearray.fromhex(nstr+"0")
	# # # print(byarray)
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
	print("=============")
#############################################

	dirs=os.listdir("./")
	for file in dirs:
		if file.endswith(".tiff"):

				img=io.imread(file)
				try:
					infosstr=endstr(img,png_wid)
					infolist=infosstr.split(",")
					if infolist[1]==	m_ts	:
						# print(file)
						os.rename(file,m_ts+"="+infolist[0]+"="+".tiff")
				except:
					oooo=1


	
	
	# # for kk in range(1,numk+1):
		# # img=io.imread( )
		
		# # img = io.imread(pngname)
	# # # # # print(img)
	
	# # bb=''
	# # for j in range(infosize):
 
		# # addnew=hex(int(img[0,j]))
		# # bb=bb+addnew




	# # nstr=bb.replace("0x","").rstrip("0")
	# # # # # print(nstr)
	# # try:
		# # byarray=bytearray.fromhex(nstr)
	# # except:
		# # byarray=bytearray.fromhex(nstr+"0")
	# # # # # print(byarray)
	# # infostr=byarray.decode("utf8")             
	# # infolist=infostr.split(",")
	
	
###############################################
	
	
	
	
	
	# # # print(len(infolist))
	
	bkname="bk_"+filename
	f=open(bkname,"wb")
	
	
	
	i=0
	j=0
	# k=1


	# # # # k<numk
	count=0
	for k in range(1,numk):
		print(k)
		img = io.imread(m_ts+"="+str(k)+"="+".tiff")
		for i in range( png_wid):
			for j in range( png_wid):
				

				
				# # print("i+\n")
				# # print(img[i,j])

				tmpbyte2= math.floor(  int(img[i,j])/2**8 )  
				tmpbyte1= math.floor(  int(img[i,j])%2**8 )  

				# print(hex(tmpbyte1))
				# print(hex(tmpbyte2))
				
				# print(hex(img[i,j]))

				# time.sleep(1)

				mywrite(f,tmpbyte1)
				mywrite(f,tmpbyte2)
	
	n2num1=num1-(numk-1)*png_wid*png_wid

	numi=math.floor(   n2num1/png_wid  )
	numj=math.floor(   n2num1%png_wid  )
	
	# # # print("(*****")
	# # # print(numi)
	# # # print(numj)
	# # # # # # k=numk
	for k in range(numk,numk+1):
		img = io.imread(m_ts+"="+str(k)+"="+".tiff")
		for i in range( numi):
			for j in range( png_wid):
				
				# # # print("\n[")
				# # # print(img[i,j])
				# # # print("]")
				
				
				# # # print(int(img[i,j])/2**8)
			
				tmpbyte2= math.floor(  int(img[i,j])/2**8 )  
				tmpbyte1= math.floor(  int(img[i,j])%2**8 )  
				
				# print(tmpbyte1)
				# print(tmpbyte2) 
				

				
				mywrite(f,tmpbyte1)
				mywrite(f,tmpbyte2)
				
		for i in range( numi,numi+1):
			for j in range( numj):
				
				# print("ii+++\n")
				# print(img[i,j])
				tmpbyte2= math.floor(  int(img[i,j])/2**8 )  
				tmpbyte1= math.floor(  int(img[i,j])%2**8 )  

				mywrite(f,tmpbyte1)
				mywrite(f,tmpbyte2)


		
		if num2>0:
			if k<numk:
				
				img = io.imread(m_ts+"="+str(numk)+"="+".tiff")
				
				tmpbyte=hex(img[0,0])
				mywrite(f,tmpbyte)

			else:


				tmpbyte=math.floor(  int(img[numi,numj])%2**8 )  
				mywrite(f,tmpbyte)


 

	f.close()	
		
	
	
	
    
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
	print(e)

if sys.argv[2]=="recover":
	print("ok")
	png2file(sys.argv[1])
				
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
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

