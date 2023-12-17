from DrissionPage import WebPage
import time



urla="https://space.bilibili.com/485156027/fans/follow"


page = WebPage()
page.get(urla)

input("开始取关")

for i in range(200):
    print(i)
    page.get("https://space.bilibili.com/485156027/fans/follow")
    time.sleep(2)
    page.run_js('''$(".be-dropdown-item:contains('取消关注')").click()''')
    time.sleep(3)

