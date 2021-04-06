from mitmproxy import ctx
import requests

def request(flow): 
    # data = {
    # "url": str(flow.request.url),
    # 'method': str(flow.request.method),
    # 'get_text': str(flow.request.get_text),
    # }

    # json_data = json.dumps(data)
    # fp = open('D:/66666.json', 'a+', encoding='utf-8')
    # fp.write(json_data + '\n')
    pass

def response(flow):
    # if flow.request.host == "mp.weixin.qq.com":
	if "mp.weixin.qq.com/mp" in flow.request.url :
		data = str(flow.response.text)
		fp = open('weixin_mitm.txt', 'a+', encoding='utf-8')
		if flow.request.query.get('__biz'):
			fp.write("\nbiz:"+flow.request.query.get('__biz')+"\n\n") 
		
		fp.write(data + '\n')
		fp.close()

addons = [
	response()
]