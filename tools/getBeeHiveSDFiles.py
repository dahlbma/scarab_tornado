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
    sFileName = sSchema + '_' + sTable + '.sdf'
    print("Exporting " + sFileName)
    try:
        os.remove(sFileName)
    except:
        pass
    
    sST = sSchema + '.' + sTable
    sdfile = open(sFileName, 'a')
    sCols = ','.join(lCols)
    sSql = 'select ' + sCols + ' from ' + sST + ' where rownum < 3000000'
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

getTable('chemspec', 'jcmol_moltable', ['regno', 'molfile'])
getTable('bcpvs', 'jcmol_moltable', ['compound_id', 'molfile'])
getTable('bcpvs', 'jcsalt_moltable', ['suffix', 'molfile'])
getTable('bcpvs', 'jcsepmol_moltable', ['compound_id', 'molfile'])

quit()
getTable('screenmol', 'jcscreenmol_moltable', ['cd_id', 'dummy_id', 'cd_formula', 'cd_molweight', 'cd_timestamp', 'molfile'])
