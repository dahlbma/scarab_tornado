import requests
import json

baseUrl = 'https://esox3.scilifelab.se/vialdb/'

def chemRegAddMolFile(dict, token):
    r = requests.post(f'{baseUrl}api/chemRegAddMol',
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
        cleanList.append(' ')
    for i in res:
        cleanList.append(i[0])
    return cleanList

def getDatabase():
    r = requests.get(f'{baseUrl}getDatabase')
    res = listify(r, False)
    return res

def searchValue(target, value, token):
    r = requests.get(f'{baseUrl}api/search',
                     params={'column': target,
                             'value': value},
                     headers={'token': token})
    cleanList = listify(r, False)
    return cleanList

def updateValue(target, value, token, regno):
    r = requests.put(f'{baseUrl}api/update',
                     params={'column': target,
                             'value': value,
                             'regno': regno},
                     headers={'token': token})

def createNewRegno(regno, token):
    r = requests.put(f'{baseUrl}api/createRegno',
                     params={'regno': regno},
                     headers={'token': token})

def updateBatch(token, regno, sBatch):
    r = requests.put(f'{baseUrl}api/updateRegnoBatch',
                     params={'regno': regno, 'batch': sBatch},
                     headers={'token': token})
    if r.status_code != 200:
        return False
    else:
        return True

def deleteRegno(regno, token):
    r = requests.put(f'{baseUrl}api/deleteRegno',
                     params={'regno': regno},
                     headers={'token': token})

def getTextColumn(token, column, regno):
    r = requests.get(f'{baseUrl}api/getTextColumn',
                     params={'column': column,
                             'regno': regno},
                     headers={'token': token})
    cleanList = listify(r, False)
    return cleanList[0]
    
def getColComboData(token, column):
    r = requests.get(f'{baseUrl}api/getColComboData',
                     params={'column': column},
                     headers={'token': token})
    cleanList = listify(r)
    return cleanList

def getLibraryName(token, library_id):
    r = requests.get(f'{baseUrl}api/getLibraryName',
                     params={'library_id': library_id},
                     headers={'token': token})
    cleanList = listify(r, False)
    if cleanList == []:
        cleanList = ' '
    else:
        cleanList = cleanList[0]
    return cleanList

def createLibrary(token, library_name, supplier):
    r = requests.put(f'{baseUrl}api/createLibrary',
                     params={'library_name': library_name,
                             'supplier': supplier},
                     headers={'token': token})
    if r.status_code != 200:
        return False, r.content
    else:
        return True, r.content
    
def createSupplier(token, supplier):
    r = requests.put(f'{baseUrl}api/createSupplier',
                     params={'supplier': supplier},
                     headers={'token': token})
    if r.status_code != 200:
        return False, r.content
    else:
        return True, r.content

def getNextRegno(token):
    r = requests.get(f'{baseUrl}api/getNextRegno',
                     headers={'token': token})
    res = r.content.decode()
    return str(res)

def getSdfSequence(token):
    r = requests.get(f'{baseUrl}api/getNextSdfSequence',
                     headers={'token': token})
    res = r.content.decode()
    return res

def getMolFile(token, regno):
    r = requests.get(f'{baseUrl}api/getMolfile', 
                     params={'regno': regno}, 
                     headers={'token': token})
    res = r.content.decode()
    return res
    
def createMolImage(token, regno):
    r = requests.get(f'{baseUrl}api/createMolImage',
                     params={'regno': regno},
                     headers={'token': token})
    res = r.content.decode()
    return res

def postMolfile(token, molfile):
    r = requests.post(f'{baseUrl}api/loadMolfile',
                      headers={'token': token},
                      files=molfile)
    return r

def login(user, password, database):
    r = requests.post(f'{baseUrl}login',
                      data = {'username':user,
                              'password':password,
                              'database':database})
    return r

def getVersion():
    r = requests.get(f'{baseUrl}getVersionData') # get file version
    return r

def getChemRegBinary(os_name):
    r = requests.get(f'{baseUrl}getChemRegBin',
                     data={'os_name':os_name},
                     stream=True) # fetch chemreg dist
    return r

def getMolImage(regno):
    r = requests.get(f'{baseUrl}mols/{regno}.png')
    res = r.content
    return res

def getLastBatchOfEln(token, sEln):
    r = requests.get(f'{baseUrl}api/getLastBatchFromEln',
                     params={'eln': sEln},
                     headers={'token': token})
    res = r.content.decode()
    return int(res)

def getCanonicSmiles(token, smiles):
    r = requests.get(f'{baseUrl}api/getCanonicSmiles',
                     params={'smiles': smiles},
                     headers={'token': token})
    return r.content.decode()

def createSalt(token, smiles):
    r = requests.put(f'{baseUrl}api/createSalt',
                     params={'smiles': smiles},
                     headers={'token': token})
    if r.status_code != 200:
        return False, r.content
    else:
        return True, r.content

def getRegnosFromSdfSequence(token, iSequence):
    r = requests.get(f'{baseUrl}api/getRegnosFromSequence',
                     params={'sdfile_sequence': iSequence},
                     headers={'token': token})
    res = listify(r, False)
    return res

def bcpvsRegCompound(token, sReg):
    r = requests.put(f'{baseUrl}api/bcpvsRegCompound',
                     params={'regno': sReg},
                     headers={'token': token})
    return r.content

def uploadBinary(token, os_name, file):
    r = requests.post(f'{baseUrl}uploadBinary',
                      data = {'os_name':os_name},
                      headers = {'token':token},
                      files = {'file':file})
    if r.status_code != 200:
        return r.content.decode(), False
    else:
        return r.content.decode(), True

def getChemRegBinary(os_name):
    r = requests.get(f'{baseUrl}getChemRegBinary/{os_name}',
                     stream=True) # fetch cello dist
    return r

def uploadVersionNo(token, ver_no):
    r = requests.post(f'{baseUrl}uploadVersionNo',
                      data = {'ver_no':ver_no},
                      headers = {'token':token}) # get file version
    if r.status_code != 200:
        return r.content.decode(), False
    else:
        return r.content.decode(), True