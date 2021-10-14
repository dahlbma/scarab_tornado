import requests
import json

def listify(data):
    res = data.content.decode()
    res = json.loads(res)
    cleanList = list()
    for i in res:
        cleanList.append(i[0])
    return cleanList
    
def getSubmitters(token):
    r = requests.get('http://esox3.scilifelab.se:8082/api/getChemists',
                     headers={'token': token})
    cleanList = listify(r)
    print(cleanList)
    return cleanList

def getProjects(token):
    r = requests.get('http://esox3.scilifelab.se:8082/api/getProjects',
                     headers={'token': token})
    cleanList = listify(r)
    return cleanList

def getCompoundTypes(token):
    r = requests.get('http://esox3.scilifelab.se:8082/api/getCompoundTypes',
                     headers={'token': token})
    cleanList = listify(r)
    return cleanList

def getProductTypes(token):
    r = requests.get('http://esox3.scilifelab.se:8082/api/getProductTypes',
                     headers={'token': token})
    cleanList = listify(r)
    return cleanList

def getNextRegno(token):
    r = requests.get('http://esox3.scilifelab.se:8082/api/getNextRegno',
                     headers={'token': token})
    res = r.content.decode()
    return res
    
