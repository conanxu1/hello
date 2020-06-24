from PyPDF2 import PdfFileReader 
import os

def num_pages(pdf_path):   


    with open(pdf_path,'rb')as f:
        pdf = PdfFileReader(f) 
 
        numb = pdf.getNumPages()


    return numb
        
        
# a=get_num_pages( '1.pdf'  )

# print(a)
# print(1)


sum=0
for parent, dirNames, fileNames in os.walk('.'):   
    for name in fileNames: 
        ext = ['pdf']            
        if name.endswith(tuple(ext)):               
            print(name+'\n')
            numm=num_pages(name)
            print(numm)
            print('\n')
            
            
            
            sum+=int(numm)
            
print(sum-351)          
            # paths_list.append(os.path.join(parent, name))
      

