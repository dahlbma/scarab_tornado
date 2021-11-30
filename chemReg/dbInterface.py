import requests
import json

def chemRegAddMolFile(dict, token):
    r = requests.post('http://esox3.scilifelab.se:8082/api/chemRegAddMol',
                      data = dict,
                      headers={'token': token})
    if r.status_code != 200:
        return False, r.content
    else:
        return True, r.content

def listify(data, addBlank=True):
    res = data.content.decode()
    res = json.loads(res)
    cleanList = list()
    if addBlank:
        cleanList.append(None)
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

def updateBatch(token, regno, sBatch):
    r = requests.put('http://esox3.scilifelab.se:8082/api/updateRegnoBatch',
                     params={'regno': regno, 'batch': sBatch},
                     headers={'token': token})
    if r.status_code != 200:
        return False
    else:
        return True

def deleteRegno(regno, token):
    r = requests.put('http://esox3.scilifelab.se:8082/api/deleteRegno',
                     params={'regno': regno},
                     headers={'token': token})

def getTextColumn(token, column, regno):
    r = requests.get('http://esox3.scilifelab.se:8082/api/getTextColumn',
                     params={'column': column, 'regno': regno},
                     headers={'token': token})
    cleanList = listify(r, False)
    return cleanList[0]
    
def getColComboData(token, column):
    r = requests.get('http://esox3.scilifelab.se:8082/api/getColComboData',
                     params={'column': column},
                     headers={'token': token})
    cleanList = listify(r)
    return cleanList

def getLibraryName(token, library_id):
    r = requests.get('http://esox3.scilifelab.se:8082/api/getLibraryName',
                     params={'library_id': library_id},
                     headers={'token': token})
    cleanList = listify(r, False)
    if cleanList == []:
        cleanList = ' '
    else:
        cleanList = cleanList[0]
    return cleanList

def createLibrary(token, library_name, supplier):
    r = requests.put('http://esox3.scilifelab.se:8082/api/createLibrary',
                     params={'library_name': library_name, 'supplier': supplier},
                     headers={'token': token})
    if r.status_code != 200:
        return False, r.content
    else:
        return True, r.content
    
def createSupplier(token, supplier):
    r = requests.put('http://esox3.scilifelab.se:8082/api/createSupplier',
                     params={'supplier': supplier},
                     headers={'token': token})
    if r.status_code != 200:
        return False, r.content
    else:
        return True, r.content

def getNextRegno(token):
    r = requests.get('http://esox3.scilifelab.se:8082/api/getNextRegno',
                     headers={'token': token})
    res = r.content.decode()
    return res

def getSdfSequence(token):
    r = requests.get('http://esox3.scilifelab.se:8082/api/getNextSdfSequence',
                     headers={'token': token})
    res = r.content.decode()
    return res

def getMolFile(token, regno):
    r = requests.get('http://esox3.scilifelab.se:8082/api/getMolfile', 
                     params={'regno': regno}, 
                     headers={'token': token})
    res = r.content.decode()
    return res
    
def createMolImage(token, regno):
    r = requests.get('http://esox3.scilifelab.se:8082/api/createMolImage',
                     params={'regno': regno},
                     headers={'token': token})
    res = r.content.decode()
    return res

def getLastBatchOfEln(token, sEln):
    r = requests.get('http://esox3.scilifelab.se:8082/api/getLastBatchFromEln',
                     params={'eln': sEln},
                     headers={'token': token})
    res = r.content.decode()
    return int(res)

def getCanonicSmiles(token, smiles):
    r = requests.get('http://esox3.scilifelab.se:8082/api/getCanonicSmiles',
                     params={'smiles': smiles},
                     headers={'token': token})
    return r.content.decode()

def createSalt(token, smiles):
    r = requests.put('http://esox3.scilifelab.se:8082/api/createSalt',
                     params={'smiles': smiles},
                     headers={'token': token})
    if r.status_code != 200:
        return False, r.content
    else:
        return True, r.content

def getRegnosFromSdfSequence(token, iSequence):
    r = requests.get('http://esox3.scilifelab.se:8082/api/getRegnosFromSequence',
                     params={'sdfile_sequence': iSequence},
                     headers={'token': token})
    res = listify(r, False)
    return res

def bcpvsRegCompound(token, sReg):
    r = requests.put('http://esox3.scilifelab.se:8082/api/bcpvsRegCompound',
                     params={'regno': sReg},
                     headers={'token': token})
    return r.content
