import pymysql



def myhand(sttr,cur):
	sql=sttr
	cur.execute(sql)
	items = cur.fetchall()
	print(items)


f=open("name_id.txt","r",encoding="utf-8")
w=f.readlines()
f.close()




MYSQL_CONFIG = {
    'host': 'localhost',  # IP地址
    'port': 3306,  # 端口
    'user': 'xuwenhan',  # 用户名
    'passwd': 'xuwenhan',  # 密码
    'db': 'bili',  # 数据库
    'charset': 'utf8',  # 编码
}

conn = pymysql.connect(**MYSQL_CONFIG)  # 数据库连接
cur = conn.cursor(pymysql.cursors.DictCursor)  # 游标对象


 
 
myhand("select * from id_name",cur)
 
 
 
 
 
# # sql = '''INSERT INTO id_name VALUE(   "'''+ans[0]+'''"   ,  "'''+ans[1]+      '''"   )'''
# # cur.execute(sql)
# # items = cur.fetchall()

for ee in w:
	ans=ee.replace("\n","").split(";")
	# print(ans)
	
	# heresql = '''insert into id_name values(   "'''+ans[0]+'''"   ,  "'''+ans[1]+      '''"   )'''
	# heresql =''' insert into id_name values(   "a","b");'''
	# print(heresql)
	# myhand(heresql,cur)
	# print(ans[0])
	# print(ans[1])
	sql = "INSERT INTO id_name VALUES (%s, %s)"
	cur.execute(sql , (ans[0],ans[1]))
	conn.commit()



	# # cur.execute(sql)


# # items = cur.fetchall()




# # print(items)
cur.close()
conn.close()








# sql = '''INSERT INTO id_name VALUE(   "'''+ "A" +'''"   ,  "'''+  "b"  +      '''"   )'''
# SQL="SELECT * FROM id_name"
# cur.execute(sql)
# items = cur.fetchall()
# print(items)




# # #创建pythonBD数据库
# # cursor.execute('CREATE DATABASE IF NOT EXISTS pythonDB DEFAULT CHARSET utf8 COLLATE utf8_general_ci;')