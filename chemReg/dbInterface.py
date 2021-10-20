import requests
import json

def listify(data, addBlank=True):
    res = data.content.decode()
    res = json.loads(res)
    cleanList = list()
    if addBlank:
        cleanList.append(' ')
    for i in res:
        cleanList.append(i[0])
    return cleanList

def searchValue(target, value, token):
    r = requests.get('http://esox3.scilifelab.se:8082/api/search',
                     params={'column': target, 'value': value},
                     headers={'token': token})
    cleanList = listify(r, False)
    return cleanList

def updateValue(target, value, token, regno):
    r = requests.put('http://esox3.scilifelab.se:8082/api/update',
                     params={'column': target, 'value': value, 'regno': regno},
                     headers={'token': token})

def createNewRegno(regno, token):
    r = requests.put('http://esox3.scilifelab.se:8082/api/createRegno',
                     params={'regno': regno},
                     headers={'token': token})
    
def deleteRegno(regno, token):
    r = requests.put('http://esox3.scilifelab.se:8082/api/deleteRegno',
                     params={'regno': regno},
                     headers={'token': token})

def getSubmitters(token):
    r = requests.get('http://esox3.scilifelab.se:8082/api/getChemists',
                     headers={'token': token})
    cleanList = listify(r)
    return cleanList

def getTextColumn(token, column, regno):
    r = requests.get('http://esox3.scilifelab.se:8082/api/getTextColumn',
                     params={'column': column, 'regno': regno},
                     headers={'token': token})
    cleanList = listify(r, False)
    return cleanList[0]
    
def getColComboData(token, regno, column):
    r = requests.get('http://esox3.scilifelab.se:8082/api/getColComboData',
                     params={'column': column},
                     headers={'token': token})
    cleanList = listify(r)
    r = requests.get('http://esox3.scilifelab.se:8082/api/getRegnoData',
                     params={'regno': regno, 'column': column},
                     headers={'token': token})
    currentValue = listify(r, False)[0]
    if currentValue != None:
        cleanList.remove(currentValue)
        cleanList.insert(0, currentValue)
    return cleanList


def getProjects(token, regno):
    r = requests.get('http://esox3.scilifelab.se:8082/api/getProjects',
                     headers={'token': token})
    cleanList = listify(r)
    r = requests.get('http://esox3.scilifelab.se:8082/api/getRegnoData',
                     params={'regno': regno, 'column': 'project'},
                     headers={'token': token})
    currentValue = listify(r, False)[0]
    if currentValue != None:
        cleanList.remove(currentValue)
        cleanList.insert(0, currentValue)
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
    
