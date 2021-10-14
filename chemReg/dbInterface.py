import requests
import json


def getSubmitters(token):
    r = requests.get('http://esox3.scilifelab.se:8082/api/getChemists',
                     headers={'token': token})
    res = r.content.decode()
    res = json.loads(res)
    cleanList = list()
    for i in res:
        cleanList.append(i[0])
    return cleanList

def getProjects(token):
    r = requests.get('http://esox3.scilifelab.se:8082/api/getProjects',
                     headers={'token': token})
    res = r.content.decode()
    res = json.loads(res)
    cleanList = list()
    for i in res:
        cleanList.append(i[0])
    return cleanList

def getCompoundTypes(token):
    r = requests.get('http://esox3.scilifelab.se:8082/api/getCompoundTypes',
                     headers={'token': token})
    res = r.content.decode()
    res = json.loads(res)
    cleanList = list()
    for i in res:
        cleanList.append(i[0])
    return cleanList

def getProductTypes(token):
    r = requests.get('http://esox3.scilifelab.se:8082/api/getProductTypes',
                     headers={'token': token})
    res = r.content.decode()
    res = json.loads(res)
    cleanList = list()
    for i in res:
        cleanList.append(i[0])
    return cleanList
