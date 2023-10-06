import requests
import json
from operator import itemgetter

# https://ifconfig.me/ip


ws=(requests.get("http://127.0.0.1:8080/get?type=HTTPS&count=all").text) 



 
json1 = json.loads(ws)

print(json1)


jsso = sorted(json1, key=itemgetter('Speed'))


for ee in jsso:

    print(ee["Ip"]\
        +" "+ee["Port"]\
        +" "+ee["Country"]\
        +" "+ee["Speed"]\
        +" "+ee["Anonymity"]\
        +" "+ee["Isp"]\
        +" "+ee["City"])
        
    print("\n\n")