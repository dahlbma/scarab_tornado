import tornado.gen
import json
import logging
import datetime
import time
import os, random, string
import re
import util
import mydb
import config
import zipfile



def getMolfile(cur, sCmpId):
    sSql = f"""select mol from bcpvs.JCMOL_MOLTABLE where compound_id = '{sCmpId}'"""
    cur.execute(sSql)
    res = cur.fetchall()
    sMol = ''
    if len(res) == 1:
        sMol = f'''{res[0][0]}
> <CIDX>
{sCmpId}

$$$$
'''
    return sMol

'''
Sample data:
EXP-14-BE7698
Helleday_SAMHD1
AA0763244 AA1963128 AA0340383

Example of ACTIVITY.tsv file:
RIDX	Pathogen_Box_Bloggs	Pathogen_Box_Bloggs	Pathogen_Box_Bloggs
CRIDX	Pathogen_Box_Bloggs	Pathogen_Box_Bloggs	Pathogen_Box_Bloggs
CRIDX_DOCID
CRIDX_CHEMBLID
CIDX	MMV161996	MMV202458	MMV676395
SRC_ID_CIDX
AIDX	1	1	1
TYPE	Inhibition	Inhibition	Inhibition
ACTION_TYPE	ANTAGONIST	ANTAGONIST	ANTAGONIST
TEXT_VALUE
RELATION	=	=	=
VALUE	0	1	61
UPPER_VALUE
UNITS	%	%	%
SD_PLUS
SD_MINUS
ACTIVITY_COMMENT	Not active	Not active	Active
ACT_ID	PB_FECH_MMV161996	PB_FECH_MMV202458	PB_FECH_MMV676395
TEOID



'''


def exportFromBatches(cur, saBatches, sRIDX, compound_record_file, molfile_file):
    for batch in saBatches:
        sMol = ''
        sSql = f"""select compound_id, "{sRIDX}", notebook_ref, compound_id from bcpvs.batch where notebook_ref = '{batch}'"""
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) == 1:
            sMol = getMolfile(res[0][0])
            compound_record_file.write('\t'.join(map(str, res[0])) + '\n')
            molfile_file.write(sMol)



def create_ACTIVITY_tsv(tRes,
                        sProject,
                        sELN,
                        sRIDX,
                        sAIDX,
                        sTARGET_TYPE,
                        sASSAY_TYPE,
                        sACTION_TYPE,
                        dir_name,
                        logg):
    '''
RIDX
CRIDX
RIDX_DOCID
CRIDX_CHEMBLID
CIDX
SRC_ID_CIDX
AIDX
TYPE
ACTION_TYPE
TEXT_VALUE
RELATION
VALUE
UPPER_VALUE
UNITS
SD_PLUS
SD_MINUS
ACTIVITY_COMMENT
ACT_ID
TEOID

tRes columns:
0    compound_id,
1    compound_batch,
2    mol,
3    inhibition,
4    activation,
5    hit_description


    '''
    iDataCol = 3
    sType = 'Inhibition'
    
    if tRes[0][3] == None and tRes[0][4] != None:
        iDataCol = 4
        sType = 'Activation'

    # sRIDX = 
    sCRIDX = sRIDX
    sCRIDX_DOCID = ''
    sCRIDX_CHEMBLID = ''
    # sCIDX = 
    sSRC_ID_CIDX = ''
    # sAIDX = 
    sTYPE = sType
    sACTION_TYPE = sACTION_TYPE
    sTEXT_VALUE = ''
    sRELATION = '='
    # sVALUE = 
    sUPPER_VALUE = ''
    sUNITS = '%'
    sSD_PLUS = ''
    sSD_MINUS = ''
    # sACTIVITY_COMMENT = 
    sACT_ID = ''
    sTEOID = sELN

        
    activity_tsv_name = dir_name + "/ACTIVITY.tsv"
    with open(activity_tsv_name, 'w') as activity_tsv_file:
        sHeader = (
            'RIDX',
            'CRIDX',
            'CRIDX_DOCID',
            'CRIDX_CHEMBLID',
            'CIDX',
            'SRC_ID_CIDX',
            'AIDX',
            'TYPE',
            'ACTION_TYPE',
            'TEXT_VALUE',
            'RELATION',
            'VALUE',
            'UPPER_VALUE',
            'UNITS',
            'SD_PLUS',
            'SD_MINUS',
            'ACTIVITY_COMMENT',
            'ACT_ID',
            'TEOID')
        activity_tsv_file.write('\t'.join(map(str, sHeader)) + '\n')

        for row in tRes:
            sCIDX = row[0] # COMPOUND_ID
            sValue = row[iDataCol]
            sActivityComment = row[5]
            activity_tsv_file.write(f'''{sRIDX}\t{sCRIDX}\t{sCRIDX_DOCID}\t{sCRIDX_CHEMBLID}\t{sCIDX}\t{sSRC_ID_CIDX}\t{sAIDX}\t{sTYPE}\t{sACTION_TYPE}\t{sTEXT_VALUE}\t{sRELATION}\t{sValue}\t{sUPPER_VALUE}\t{sUNITS}\t{sSD_PLUS}\t{sSD_MINUS}\t{sActivityComment}\t{sACT_ID}\t{sTEOID}\n''')

        
def exportFromElnProject(cur,
                         saCompounds,
                         sProject,
                         sELN,
                         sRIDX,
                         sAIDX,
                         sTARGET_TYPE,
                         sASSAY_TYPE,
                         sACTION_TYPE,
                         compound_record_file,
                         molfile_file,
                         dir_name,
                         logg):
    sComp = ''
    sNotTheseCompounds = ''
    if len(saCompounds) == 1:
        sComp = f''' '{saCompounds[0]}' '''
    elif len(saCompounds) > 1:
        for compound in saCompounds:
            if sComp == '':
                sComp = f''' '{compound}' '''
            else:
                sComp += f''', '{compound}' '''

    if sComp != '':
        sNotTheseCompounds = f''' and a.compound_id not in ({sComp}) '''

    sSql = f'''select a.compound_id,
    compound_batch,
    bin2mol( moldepict(mol) ),
    ROUND(inhibition * 100, 4) AS inhibition,
    ROUND(activation * 100, 4) AS activation,
    CASE 
        WHEN hit = 0 THEN 'Not active'
        WHEN hit = 1 THEN 'Active'
    END AS hit_description

    from assay.lcb_sp a, bcpvs.JCMOL_MOLTABLE m
    where a.compound_id = m.compound_id
    and a.project = '{sProject}'
    and a.eln_id = '{sELN}'
    {sNotTheseCompounds}
    order by a.compound_id
    '''
    try:
        cur.execute(sSql)
    except:
        logg.error(f"{sSql}")

    res = cur.fetchall()
    sPrevCompId = ''
    for row in res:
        sCmpId = row[0]
        sBatch = row[1]
        sMol = row[2]
        sMolfile = ''
        
        if len(sMol) > 5 and sPrevCompId != sCmpId:
            sMolfile = f'''{sMol.decode('utf-8')}
> <CIDX>
{sCmpId}

$$$$
'''

            lines = sMolfile.split('\n')
            # Replace the first line with an empty string
            lines[0] = ''
            # Join the lines back into a single string
            sMolfile = '\n'.join(lines)
            
            compound_record_file.write(f'''{sCmpId}\t{sRIDX}\t{sBatch}\t{sCmpId}\n''')
            molfile_file.write(sMolfile)
        sPrevCompId = sCmpId
            
    create_ACTIVITY_tsv(res,
                        sProject,
                        sELN,
                        sRIDX,
                        sAIDX,
                        sTARGET_TYPE,
                        sASSAY_TYPE,
                        sACTION_TYPE,
                        dir_name,
                        logg)


