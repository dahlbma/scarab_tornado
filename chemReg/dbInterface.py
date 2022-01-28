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
        cleanList.append(None)
    for i in res:
        cleanList.append(i[0])
    return cleanList

def getDatabase():
    r = requests.get(f'{baseUrl}getDatabase')
    res = listify(r, False)
    return res

def searchValue(target, value, token, database):
    r = requests.get(f'{baseUrl}api/search',
                     params={'column': target,
                             'value': value,
                             'database': database},
                     headers={'token': token})
    cleanList = listify(r, False)
    return cleanList

def updateValue(target, value, token, regno, database):
    r = requests.put(f'{baseUrl}api/update',
                     params={'column': target,
                             'value': value,
                             'regno': regno,
                             'database': database},
                     headers={'token': token})

def createNewRegno(regno, token, database):
    r = requests.put(f'{baseUrl}api/createRegno',
                     params={'regno': regno, 'database': database},
                     headers={'token': token})

def updateBatch(token, regno, sBatch, database):
    r = requests.put(f'{baseUrl}api/updateRegnoBatch',
                     params={'regno': regno, 'batch': sBatch, 'database': database},
                     headers={'token': token})
    if r.status_code != 200:
        return False
    else:
        return True

def deleteRegno(regno, token, database):
    r = requests.put(f'{baseUrl}api/deleteRegno',
                     params={'regno': regno, 'database': database},
                     headers={'token': token})

def getTextColumn(token, column, regno, database):
    r = requests.get(f'{baseUrl}api/getTextColumn',
                     params={'column': column,
                             'regno': regno,
                             'database': database},
                     headers={'token': token})
    cleanList = listify(r, False)
    return cleanList[0]
    
def getColComboData(token, column, database):
    r = requests.get(f'{baseUrl}api/getColComboData',
                     params={'column': column, 'database': database},
                     headers={'token': token})
    cleanList = listify(r)
    return cleanList

def getLibraryName(token, library_id, database):
    r = requests.get(f'{baseUrl}api/getLibraryName',
                     params={'library_id': library_id, 'database': database},
                     headers={'token': token})
    cleanList = listify(r, False)
    if cleanList == []:
        cleanList = ' '
    else:
        cleanList = cleanList[0]
    return cleanList

def createLibrary(token, library_name, supplier, database):
    r = requests.put(f'{baseUrl}api/createLibrary',
                     params={'library_name': library_name,
                             'supplier': supplier,
                             'database': database},
                     headers={'token': token})
    if r.status_code != 200:
        return False, r.content
    else:
        return True, r.content
    
def createSupplier(token, supplier, database):
    r = requests.put(f'{baseUrl}api/createSupplier',
                     params={'supplier': supplier, 'database': database},
                     headers={'token': token})
    if r.status_code != 200:
        return False, r.content
    else:
        return True, r.content

def getNextRegno(token, database):
    r = requests.get(f'{baseUrl}api/getNextRegno',
                     params={'database': database},
                     headers={'token': token})
    res = r.content.decode()
    return str(res)

def getSdfSequence(token, database):
    r = requests.get(f'{baseUrl}api/getNextSdfSequence',
                     params={'database': database},
                     headers={'token': token})
    res = r.content.decode()
    return res

def getMolFile(token, regno, database):
    r = requests.get(f'{baseUrl}api/getMolfile', 
                     params={'regno': regno, 'database': database}, 
                     headers={'token': token})
    res = r.content.decode()
    return res
    
def createMolImage(token, regno, database):
    r = requests.get(f'{baseUrl}api/createMolImage',
                     params={'regno': regno, 'database': database},
                     headers={'token': token})
    res = r.content.decode()
    return res

def postMolfile(token, molfile, database):
    r = requests.post(f'{baseUrl}api/loadMolfile',
                      params={'database': database},
                      headers={'token': token},
                      files=molfile)
    return r

def login(user, password):
    r = requests.post(f'{baseUrl}login',
                      data = {'username':user,
                              'password':password})
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

def getLastBatchOfEln(token, sEln, database):
    r = requests.get(f'{baseUrl}api/getLastBatchFromEln',
                     params={'eln': sEln, 'database': database},
                     headers={'token': token})
    res = r.content.decode()
    return int(res)

def getCanonicSmiles(token, smiles, database):
    r = requests.get(f'{baseUrl}api/getCanonicSmiles',
                     params={'smiles': smiles, 'database': database},
                     headers={'token': token})
    return r.content.decode()

def createSalt(token, smiles, database):
    r = requests.put(f'{baseUrl}api/createSalt',
                     params={'smiles': smiles, 'database': database},
                     headers={'token': token})
    if r.status_code != 200:
        return False, r.content
    else:
        return True, r.content

def getRegnosFromSdfSequence(token, iSequence, database):
    r = requests.get(f'{baseUrl}api/getRegnosFromSequence',
                     params={'sdfile_sequence': iSequence, 'database': database},
                     headers={'token': token})
    res = listify(r, False)
    return res

def bcpvsRegCompound(token, sReg, database):
    r = requests.put(f'{baseUrl}api/bcpvsRegCompound',
                     params={'regno': sReg, 'database': database},
                     headers={'token': token})
    return r.content
