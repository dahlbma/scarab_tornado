import codecs
import config
import MySQLdb
import re


def doSearch2(sMol):
    sSql = f"""
SELECT * FROM bcpvs.`JCMOL_MOLTABLE_ukey` T1
WHERE T1.molkey = UNIQUEKEY('{sMol}')
    """
    sSql = f"""
SELECT * FROM bcpvs.`JCMOL_MOLTABLE_ukey` T1
WHERE T1.molkeyct = UNIQUEKEY('{sMol}', 'cistrans')
    """
    cur.execute(sSql)
    mols = cur.fetchall()
    return mols


def getNextMolecule(sFile):
    sMol = b""
    iCount = 0
    while True:
        iCount += 1
        if iCount > 1000:
            return ""
        line = sFile.readline()
        line = line.replace(b'\r\n', b'\n')
        line = line.replace(b"'", b"")

        # RDKit can't handle empty first line in molfile
        if iCount == 1:
            line = b'id' + line

        sMol += line
        if b'$$$$' in line:
            sMol = sMol.decode(errors='replace')
            return sMol
    return ""


def getCompoundId(sMolfile):
    x = re.match('(CBK.+)', sMolfile)
    x = re.findall('(CBK.+)', sMolfile)
    try:
        sCompoundId = x[0]
    except:
        sCompoundId = ''
        print(sMolfile)
    return sCompoundId


db_connection = MySQLdb.connect(
    host=config.database['host'],
    user=config.database['user'],
    passwd=config.database['password']
)

cur = db_connection.cursor()

sdfilename = 'sdf.sdf'
sdfilename = 'jcmol_moltable_org.sdf'
#sdfilename = 'strange_mols.sdf'
#sdfilename = 'mol.sdf'
sdfilename = 'jcmol_moltable.sdf'

f = open(sdfilename, "rb")
fErr = open("error_standardized.sdf", "w")
iMols = 0
iNothingFound = 0
iFound = 0
while True:
    sMol = getNextMolecule(f)
    sCompoundId = getCompoundId(sMol)
    if sMol == b'' or sMol == '':
        break
    iMols += 1

    if iMols % 1000 == 0:
        print(iMols)

    #if iMols % 10000 == 0:
    #    break
    #    #quit()
    #    #print(iMols)

    sSql = f"""select bin2inchi(mol2bin('{sMol}','mol'), 'key')"""
    sSql = f"""select bin2inchi(mol2bin('{sMol}','mol'))"""
    cur.execute(sSql)
    res = cur.fetchall()
    try:
        #print(res[0][0].decode())
        pass
    except:
        pass
        #print(sMol)
        #quit()
    sSql = f'''select bin2smiles(mol2bin('{sMol}')) smiles'''
    cur.execute(sSql)
    res = cur.fetchall()

    try:
        sSmiles = (res[0][0]).decode()
    except:
        pass

    mols = doSearch2(sMol)
    
    if len(mols) > 0:
        #print(len(mols))
        iFound += 1
    else:
        fErr.write(sMol)
        print(sCompoundId)
        iNothingFound += 1


print(f'Not found: {iNothingFound}')
print(f'Found: {iFound}')
print(f'Total mols: {iMols}')
