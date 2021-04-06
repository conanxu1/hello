3import json
import time

def myti(now):

	timeArray = time.localtime(now)
	#print(timeArray)
	otherStyleTime = time.strftime("%Y--%m--%d %H:%M:%S", timeArray)
	return otherStyleTime


filename=input("文件名\n")
#filename="t.txt"
 

f=open(filename,"r",encoding="utf-8")
w=f.readlines()
f.close()



for ee in w:
	try:
		flag=0
		if ee.find("&quot;")>0:
			en=ee.replace("&quot;","\"").replace("\'","\n").split("\n")

			flag=1
			print(en[1])
			print("\n\nmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm")
			ee=en[1]

		data = json.loads(ee)
# =============================================================================
# 		print("\n\n")
# 		print(type(data["general_msg_list"]))
# 		
# 		print(data["general_msg_list"])
# =============================================================================
	
		if flag==0:
			data2=json.loads(data["general_msg_list"])
            
		if flag==1:
			data2=data
		


# =============================================================================
# 		print(data2)
# =============================================================================
		
		
		f=open("picurl.txt","a",encoding="utf-8")
		fh=open(filename+".html","a",encoding="utf-8")
		fd=open(filename+".txt","a",encoding="utf-8")
		for gg in data2["list"]:
		 
			print(myti(gg["comm_msg_info"]["datetime"]))
            
			f.write(gg["app_msg_ext_info"]["cover"]+"\n")
			
			fh.write("<div><a href=\""+  gg["app_msg_ext_info"]["content_url"]+      "\">"+ gg["app_msg_ext_info"]["title"] +"    "+gg["app_msg_ext_info"]["digest"]+"</a>   </div>  </br>\n")
			fh.write("<div>"+myti(gg["comm_msg_info"]["datetime"])+ gg["app_msg_ext_info"]["content"] +gg["app_msg_ext_info"]["source_url"]+"</div></br>\n")
			fh.write("<div>    <img  width=\"200\" height=\"200\"   src="+gg["app_msg_ext_info"]["cover"]+"/>     </div>")
			
			fd.write(gg["app_msg_ext_info"]["title"]+"\n"+gg["app_msg_ext_info"]["content_url"]+"\n")

			print(gg["app_msg_ext_info"]["multi_app_msg_item_list"] )
			for ggg in  gg["app_msg_ext_info"]["multi_app_msg_item_list"]:
				f.write(ggg["cover"]+"\n")
				fh.write("<div  >   <a href=\""+  ggg["content_url"]+      "\">"+ ggg["title"] +"    "+ggg["digest"]+"</a></div></br>\n")
				fh.write("<div >"+ ggg["content"] +ggg["source_url"]+"</div></br>\n")
				fh.write("<div >     <img width=\"20\" height=\"20\"   src="+ggg["cover"]+"/>         </div>")
				fd.write(ggg["title"]+"\n"+ggg["content_url"]+"\n")

                
                            
             
			
            
			print("ok")
    
		f.close()
		fh.close()     
		fd.close()  
			
		
		
		
		
		
 
            
            
# =============================================================================
#             		print(ee["cover"])
# =============================================================================
# 			
#             
#             
#             
# 			
#             
# 			fd.write(ee["title"]+"\n"+ee["link"]+"\n")
#         
#             
# 			print(ee["title"])
#             
#             
# 			print(ee["digest"])
# 			print(ee["link"])
# 			
# 			print("\n")
#             

# =============================================================================
# =============================================================================
            
            

	except:
		pass
        #print("\n\ner!!\n")
