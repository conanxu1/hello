import pymysql
 
 
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

  


 
results=myhand("select * from id_name",MYSQL_CONFIG)
 



f=open("name_id_read.txt","w",encoding="utf-8")
for ee in results:
	print("\nid:")
	print(ee["id"])
	print("\nname:")
	print(ee["name"])
    f.write(ee["id"]+";;;"+ee["name"]+"\n")
	print("\n")



f.close()
 
print(len(results))








 