import mysql.connector
import mysql
import argparse
import config

'''
## Example usage:
python register_sdfile.py chemspec chemspec_mol regno chemspec_jcmol_moltable.sdf
##
'''
parser = argparse.ArgumentParser(description='Import a SDFile into mysql')

parser.add_argument('SchemaName', metavar='SchemaName', type=str, nargs=1,
                    help='Schema name for table')
parser.add_argument('TableName', metavar='TableName', type=str, nargs=1,
                    help='Table to create')
parser.add_argument('IdCol', metavar='IdCol', type=str, nargs=1,
                    help='Identifier for the molecule')
parser.add_argument('InputFile', metavar='InputFile', type=str, nargs=1,
                    help='Input SDfile')
args = parser.parse_args()
sSDfile = args.InputFile[0]
sSchema = args.SchemaName[0]
sTable = args.TableName[0]
sMolId = args.IdCol[0]


db_connection = mysql.connector.connect(
    host=config.database['host'],
    user=config.database['user'],
    passwd=config.database['password'],
    database=sSchema
)

cur = db_connection.cursor()

def getNextMolecule(sFile):
    sMol = "'"
    iCount = 0
    while True:
        iCount += 1
        if iCount > 1000:
            return ""
        line = sFile.readline()
        line = line.replace("'", "")
        sMol += line
        if '$$$$' in line:
            return sMol[:-1] + "'"
    return ""

def createTable(sSchema, sTable, sOtherCols):
    sSql = f"""
    create table {sSchema}.{sTable}
    ( molid INTEGER PRIMARY KEY AUTO_INCREMENT,
    mol mediumblob,
    {sOtherCols}) character set latin1 ENGINE = InnoDB
    """
    cur.execute(sSql)
    '''
    sSql = f"""
    create table {sSchema}.{sTable}_key
    ( molid INTEGER PRIMARY KEY,
    molkey char(192) BINARY) character set latin1 ENGINE = InnoDB
    """
    cur.execute(sSql)

    sSql = f"""
    create table {sSchema}.{sTable}_keysim
    ( molid INTEGER PRIMARY KEY,
    molkey char(192) BINARY) character set latin1 ENGINE = InnoDB
    """
    cur.execute(sSql)

    sSql = f"""
    create table {sSchema}.{sTable}_ukey
    ( molid INTEGER PRIMARY KEY,
    molkey mediumblob) character set latin1 ENGINE = InnoDB
    """
    cur.execute(sSql)
    '''
    
def dropTables(sSchema, sTable):
    sSql = f"""
    drop table {sSchema}.{sTable}
    """
    try:
        cur.execute(sSql)
    except:
        print('Cant drop table')
        pass

    sSql = f"""
    drop table {sSchema}.{sTable}_key
    """
    try:
        cur.execute(sSql)
    except:
        print('no _key table')

    sSql = f"""
    drop table {sSchema}.{sTable}_keysim
    """
    try:
        cur.execute(sSql)
    except:
        print('no _keysim table')

    sSql = f"""
    drop table {sSchema}.{sTable}_ukey
    """
    try:
        cur.execute(sSql)
    except:
        print('no _ukey table')
        
def insertMolecule(sSchema, sTable, sMol, sCols, sValues):
    sSql = f"""
    insert into {sSchema}.{sTable}
    (mol, {sCols}) values (mol2bin({sMol}), {sValues})
    """
    cur.execute(sSql)
    db_connection.commit()
       
def getValues(sMol):
    sVals = ""
    saMols = sMol.splitlines()
    for i in range(len(saMols)):
        if '> <' in saMols[i]:
            sVals += "'" + saMols[i+1] + "', "

    return sVals[:-2]

file = open(sSDfile,'r')

# Extra columns
sCols = f"""{sMolId} varchar(20),
            cd_molweight double,
            cd_timestamp datetime"""

dropTables(sSchema, sTable)
createTable(sSchema, sTable, sCols)

sCols = f"""{sMolId}, cd_molweight, cd_timestamp"""

iCount = 0
while True:
    sMol = getNextMolecule(file)
    if sMol == "":
        # All molecules are done
        break
    if len(sMol) > 10:
        sValues = getValues(sMol)
        insertMolecule(sSchema, sTable, sMol, sCols, sValues)
    iCount += 1

quit()
sSql = f"""
    insert into {sSchema}.{sTable}_key
    select molid, fp(mol, 'sss') as molkey from {sSchema}.{sTable}
"""
cur.execute(sSql)

sSql = f"""
    insert into {sSchema}.{sTable}_keysim
    select molid, fp(mol, 'sim') as molkey from {sSchema}.{sTable}
"""
cur.execute(sSql)

sSql = f"""
    insert into {sSchema}.{sTable}_ukey
    select molid, uniquekey(mol) as molkey from {sSchema}.{sTable}
"""
cur.execute(sSql)
db_connection.commit()

