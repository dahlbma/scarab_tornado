# -*- coding: utf-8 -*-
import ast
import datetime
import time

import os, random, string
import csv
import re
import cx_Oracle
import codecs
import mysql
import mysql.connector
import config

mysql_connection = mysql.connector.connect(
    host=config.database['host'],
    user=config.database['user'],
    passwd=config.database['password']
)


mysql_cur = mysql_connection.cursor()

os.environ['NLS_LANG'] ='AMERICAN_AMERICA.WE8ISO8859P1'
dsn_tns = cx_Oracle.makedsn(config.oracleHost, config.oraclePort, config.oracleInstance, 'latin1db')

con = cx_Oracle.connect(config.Ora['user'],
                        config.Ora['password'],
                        config.Ora['host'])

cur = con.cursor()

def getTable(sSchema, sTable, lCols):
    sFileName = f'1{sSchema}_{sTable}{sTable}.sdf'
    print("Exporting " + sFileName)
    try:
        os.remove(sFileName)
    except:
        pass
    
    sST = sSchema + '.' + sTable
    sdfile = open(sFileName, 'a')
    sCols = ','.join(lCols)
    sSql = f'''select {sCols} from {sST} where compound_id in (
'CBK664927',
'CBK664926',
'CBK664925')
 '''
    cur.execute(sSql)
    s = True
    while s:
        try:
            ll = cur.fetchone()
            sMol = ll[len(ll)-1].read()
            iCols = len(ll)
            if len(sMol) > 10:
                sdfile.write(sMol + '\n')
                for i in range(iCols-1):
                    sdfile.write('> <' + lCols[i] + '>\n')
                    sdfile.write(str(ll[i]) + '\n')
                sdfile.write('\n$$$$\n')

        except Exception as e:
            if "has no len" not in str(e):
                print(str(e))
            s = False
            break

getTable('bcpvs', 'jcmol_moltable', ['compound_id', 'molfile'])
getTable('bcpvs', 'jcsepmol_moltable', ['compound_id', 'molfile'])




'''
'CBK000212D',
'CBK000953',
'CBK000469C',
'CBK000486C',
'CBK000490C',
'CBK000491C',
'''
