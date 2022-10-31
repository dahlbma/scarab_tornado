import requests
import json

baseUrl = 'https://esox3.scilifelab.se/chemreg/'

def listify(data, addBlank=True):
    res = data.content.decode()
    res = json.loads(res)
    cleanList = list()
    if addBlank:
        cleanList.append(' ')
    for i in res:
        cleanList.append(i[0])
    if cleanList == []:
        cleanList = [None]
    return cleanList

def getDatabase():
    r = requests.get(f'{baseUrl}getDatabase')
    res = listify(r, False)
    return res

def getStructure(token, compound_id):
    r = requests.get(f'{baseUrl}api/getNextSdfSequence',
                     params={'compound_id': compound_id},
                     headers={'token': token})
    res = r.content.decode()
    return res

def getMolFile(token, regno):
    r = requests.get(f'{baseUrl}api/getMolfile', 
                     params={'regno': regno},
                     headers={'token': token})
    res = r.content.decode()
    return res


def getForwardRegno(token, regno):
    r = requests.get(f'{baseUrl}api/getForwardRegno', 
                     params={'regno': regno},
                     headers={'token': token})
    res = r.content.decode()
    return res


def getBackwardRegno(token, regno):
    r = requests.get(f'{baseUrl}api/getBackwardRegno',
                     params={'regno': regno},
                     headers={'token': token})
    res = r.content.decode()
    return res


def getRegnoFromCompound(token, sCmpId):
    r = requests.get(f'{baseUrl}api/getRegnoFromCompound', 
                     params={'compound_id': sCmpId},
                     headers={'token': token})
    res = r.content.decode()
    return res

def getCompoundFromRegno(token, sRegno):
    r = requests.get(f'{baseUrl}api/getCompoundFromRegno', 
                     params={'regno': sRegno},
                     headers={'token': token})
    res = r.content.decode()
    res = res.replace('"', '')
    return res

def getMolFileBcpvs(token, compound_id):
    r = requests.get(f'{baseUrl}api/getMolfileBcpvs', 
                     params={'compound_id': compound_id},
                     headers={'token': token})
    res = r.content.decode()
    return res

def createBcpvsMolImage(token, compound_id):
    r = requests.get(f'{baseUrl}api/createBcpvsMolImage',
                     params={'compound_id': compound_id},
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

def createMolImageFromMolfile(token, molfile):
    dValues = {
        "mol_id": 'new_structure',
        "molfile": molfile
    }
    r = requests.post(f'{baseUrl}api/createMolImageFromMolfile',
                      data = dValues,
                      headers={'token': token})
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

def getMolImage(regno):
    r = requests.get(f'{baseUrl}mols/{regno}.png')
    res = r.content
    return res

def getCanonicSmiles(token, smiles):
    r = requests.get(f'{baseUrl}api/getCanonicSmiles',
                     params={'smiles': smiles},
                     headers={'token': token})
    return r.content.decode()
