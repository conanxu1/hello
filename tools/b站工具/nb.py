import dryscrape
from bs4 import BeautifulSoup
import pymysql
import time
import webkit_server


def duplicate_removal(datas, condition, model="key"):
    """
    :param datas: 准备去重的数据，格式[{},{}...]
    :param condition: 去重参考的键值，需要数据里面有这些key
    :param model: 去重模式，key模式为去重参考的key值；notkey模式为去重不参考的key值。相反关系。
    :return: 去重后的数据，格式[{},{}...]
    """
    def flags(keys, data):
        tmp_dic = {}
        for key in keys:
            tmp_dic.update({key: data.get(key)})
        return tmp_dic

    removal_data = []
    values = []
    if datas:
        if model == "key":
            keys = condition
        elif model == "notkey":
            keys = [key for key in datas[0].keys() if key not in condition]
        else:
            raise ValueError("传入的model值错误，无法匹配")
        for data in datas:
            if flags(keys, data) not in values:
                removal_data.append(data)
                values.append(flags(keys, data))

    return removal_data

MYSQL_CONFIG = {
    'host': 'localhost',  # IP地址
    'port': 3306,  # 端口
    'user': 'xuwenhan',  # 用户名
    'passwd': 'xuwenhan',  # 密码
    'db': 'bili',  # 数据库
    'charset': 'utf8',  # 编码
}
 
 
def myhand(sttr,mysql_config):
	conn = pymysql.connect(**mysql_config)  # 数据库连接
	cur = conn.cursor(pymysql.cursors.DictCursor)  # 游标对象
	sql=sttr
	cur.execute(sql)
	items = cur.fetchall()
	cur.close()
	conn.close()
	return items
    
    
def myupdate(sttr,mysql_config,data):
    conn = pymysql.connect(**mysql_config)  # 数据库连接
    cur = conn.cursor(pymysql.cursors.DictCursor)  # 游标对象
    sql=sttr
    cur.execute(sql,data)
    conn.commit()
    cur.close()
    conn.close()

    
    

results=myhand("select * from copy_id_name",MYSQL_CONFIG)

key = ["id"]
print(len( duplicate_removal(results, key, model="key")))





def get_rep(myurl):
    
    
    
    global session
    session.set_timeout(40)
    session.set_attribute('auto_load_images',False)
    session.clear_cookies()
    
    
    
    try:
        session.visit(myurl) 
        time.sleep(10)
        response =session.body() 
    except:
        response="<b></b>"
    session.reset()
    
    
    return response



# UPDATE  copy_id_name
# SET up=0

dryscrape.start_xvfb()
session = dryscrape.Session()
 
for ee in results:
    
    print("\nid:")
    print(ee["id"])
    print("\nname:")
    print(ee["name"])
    
    my_url = "https://space.bilibili.com/"+str(ee["id"])+"/video"
    
    
    
    print(ee["up"])
    if ee["up"]==0:
    
        

        
        response=get_rep(my_url)
        
        for ii in range(2):
            rp=get_rep(my_url)
            response=rp+response
        print(response)    
        
        soup = BeautifulSoup(response,features="lxml")
        # list1=soup.find_all("ul",class_="list-list")
        
        list1=soup.find_all("li",class_="new")
       
        
        list1=set(list1)
        list1=list(list1)
        print(list1)
        
        f=open("data.html","a+",encoding="utf-8")
        
        if list1 != []:
            f.write('''<a href="http://space.bilibili.com/'''+ee["id"]+'''">'''+ee["name"]+"</a>")
            f.write('''<hr style="height:8px;border:none;border-top:3px double #9492ff;"/>''')
        
        
        for ff in list1:
            f.write('''<ul  style="column-count:2;column-width: 100px;font-size:180%;break-inside: auto!important;column-gap: 40px;" >''')
            f.write(str(ff))
            f.write("</ul>")
        f.close()


           # print(ff.prettify())
    
    myupdate("UPDATE copy_id_name SET up = 1 WHERE id = %s;",MYSQL_CONFIG,(str(ee["id"])))
    
  

        
    
    
  
	
 
print(len(results))








































