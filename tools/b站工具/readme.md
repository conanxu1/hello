### bili_pic_vid.py

获取图片或者视频地址xz

---

### bili_space.py

获取用户以及地址

---

### bili_quguan.py

批量取关脚本

---

### mysql_handle.py

上传mysql





union去重创建新表

```sql
create table `id_name20231217`
as
(select * from `id_name`)
union 
(select * from `id_name2`)
ORDER BY `id` DESC

```