import requests
import json

baseUrl = 'https://esox3.scilifelab.se/chemreg/'

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
    if cleanList == []:
        cleanList = [None]
    return cleanList

def getDatabase():
    r = requests.get(f'{baseUrl}getDatabase')
    res = listify(r, False)
    return res

def createNewRegno(regno, token):
    r = requests.put(f'{baseUrl}api/createRegno',
                     params={'regno': regno},
                     headers={'token': token})


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

def getMolFileBcpvs(token, compound_id):
    r = requests.get(f'{baseUrl}api/getMolfileBcpvs', 
                     params={'compound_id': compound_id},
                     headers={'token': token})
    res = r.content.decode()
    print(res)
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

def getMolImage(regno):
    r = requests.get(f'{baseUrl}mols/{regno}.png')
    res = r.content
    return res

def getCanonicSmiles(token, smiles):
    r = requests.get(f'{baseUrl}api/getCanonicSmiles',
                     params={'smiles': smiles},
                     headers={'token': token})
    return r.content.decode()
