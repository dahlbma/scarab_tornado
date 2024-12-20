import tornado.web
import json
import MySQLdb
from auth import jwtauth
import config
import logging
import rdkit.Chem
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from molmass import Formula
from io import StringIO
import sys
import codecs
import re
import os
import mydb
import chembl_export
import zipfile
import random

logger = logging.getLogger(__name__)

db = mydb.disconnectSafeConnect()
cur = db.cursor()


def getDatabase(parent):
    data = parent.request.headers['Token']
    jsonData = json.loads(data)
    database = jsonData['database']
    if database == 'Live':
        return 'chem_reg', 'bcpvs'
    else:
        return 'chem_reg_test', 'bcpvs_test'

def res2json():
    result = [list(i) for i in cur.fetchall()]
    return json.dumps(result)


def checkUniqueStructure(stdSMILES, bcpvsDB):
    sSql = f"""
SELECT compound_id FROM {bcpvsDB}.`compound` T1
WHERE smiles_std = '{stdSMILES}'
    """
    cur.execute(sSql)
    mols = cur.fetchall()
    return mols


def getNewRegno(chemregDB):
    sSql = f"select pkey from {chemregDB}.regno_sequence"
    cur.execute(sSql)
    pkey = cur.fetchall()[0][0] +1
    sSql = f"update {chemregDB}.regno_sequence set pkey={pkey}"
    cur.execute(sSql)
    return pkey

def getNewSaltNumber(chemregDB):
    sSql = f"select pkey from {chemregDB}.salts_sequence"
    cur.execute(sSql)
    pkey = cur.fetchall()[0][0] +1
    sSql = f"update {chemregDB}.salts_sequence set pkey={pkey}"
    cur.execute(sSql)
    return pkey

def getNewCompoundId(bcpvsDB):
    sSql = f"select pkey from {bcpvsDB}.compound_id_sequence"
    cur.execute(sSql)
    pkey = cur.fetchall()[0][0] +1
    sSql = f"update {bcpvsDB}.compound_id_sequence set pkey={pkey}"
    cur.execute(sSql)
    return pkey

def getAtomicComposition(saComp):
    sComp = f""
    for atom in saComp:
        sComp += f'{atom[0]} {round(atom[3] * 100, 2)}% '
    return sComp

def createPngFromMolfile(regno, molfile):
    m = Chem.MolFromMolBlock(molfile, removeHs=False, sanitize=True)
    try:
        Draw.MolToFile(m, f'mols/{regno}.png', kekulize=True, size=(280, 280))
    except:
        logger.error(f"regno {regno} is nostruct")


def isItNewStructure(self, molfile):
    chemregDB, bcpvsDB = getDatabase(self)
    molfile, stdSMILES = cleanStructureRDKit(molfile)
    if molfile == False:
        return False

    mols = checkUniqueStructure(stdSMILES, bcpvsDB)

    if mols == False:
        return False
    sStatus = 'oldMolecule'
    if len(mols) == 0:
        sStatus = 'newMolecule'
        compound_id = ''
    else:
        compound_id = mols[0][0]
    return compound_id


def cleanStructureRDKit(molfile):
    # Function to clean molecule from charges etc.
    mol = Chem.MolFromMolBlock(molfile, removeHs=False, sanitize=True)
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = False
    params.tautomerRemoveBondStereo = False
    params.tautomerRemoveIsotopicHs = False
    try:
        clean_mol = rdMolStandardize.Cleanup(mol, params)
    except Exception as e:
        logger.error(f"Failed to standardize {str(e)}")
        return False, False
    try:
        parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol, params)
    except:
        logger.error('Failed with rdMolStandardize.FragmentParent, skipping molecule')
        return False, False

    uncharger = rdMolStandardize.Uncharger()
    uncharged_parent_clean_mol = uncharger.uncharge(parent_clean_mol)

    te = rdMolStandardize.TautomerEnumerator(params) # idem
    taut_uncharged_parent_clean_mol = te.Canonicalize(uncharged_parent_clean_mol)
    mV4 = Chem.MolToSmiles(taut_uncharged_parent_clean_mol)
    mV4 = mV4.replace('\\', '\\\\')

    sSql = f'''select CONVERT(bin2mol(mol2bin('{mV4}' , 'smiles')) USING UTF8)'''
    cur.execute(sSql)
    saRes = cur.fetchall()

    if len(saRes) == 1:
        return saRes[0][0], mV4
    else:
        return False, False

def addStructure(database, molfile, newRegno, idColumnName):
    #####
    # Add the molecule to the structure tables
    sSql = f"""
    insert into {database} (mol, {idColumnName})
    value
    ('{molfile}', '{newRegno}')
    """
    cur.execute(sSql)

    sSql = f"""
    insert into {database}_MOL (mol, {idColumnName})
    value
    (mol2bin('{molfile}', 'mol'), '{newRegno}')
    """
    cur.execute(sSql)

    sSql = f"""
    insert into {database}_ukey select {idColumnName}, uniquekey(mol) as molkey,
    uniquekey(`mol`,'nostereo') as molkeyns,
    uniquekey(`mol`,'cistrans') as molkeyct
    from {database} where {idColumnName} = '{newRegno}'
    """
    cur.execute(sSql)

    sSql = f"""
    insert into {database}_MOL_keysim select {idColumnName}, fp(mol, 'sim') as molkey
    from {database}_MOL where {idColumnName} = '{newRegno}'
    """
    cur.execute(sSql)

    sSql = f"""
    insert into {database}_MOL_key select {idColumnName}, fp(mol, 'sss') as molkey
    from {database}_MOL where {idColumnName} = '{newRegno}'
    """
    cur.execute(sSql)
    #
    #####

def registerNewCompound(bcpvsDB,
                        compound_id_numeric,
                        molfile,
                        stdSmiles,
                        mf,
                        sep_mol_monoiso_mass,
                        ip_rights = '',
                        compound_name = ''):
    compound_id = f'CBK{compound_id_numeric}'
    sSql = f'''
    insert into {bcpvsDB}.compound (
    compound_id,
    compound_id_numeric,
    created_date,
    mf,
    ip_rights,
    sep_mol_monoiso_mass,
    smiles_std)
    values (
    '{compound_id}',
    {compound_id_numeric},
    now(),
    '{mf}',
    '{ip_rights}',
    {sep_mol_monoiso_mass},
    '{stdSmiles}')
    '''
    cur.execute(sSql)
    return compound_id

def registerNewBatch(bcpvsDB,
                     compound_id,
                     chemreg_regno,
                     notebook_ref,
                     suffix,
                     submitter,
                     project,
                     supplier,
                     biological_mw,
                     library_id,
                     compound_type,
                     product_type,
                     supplier_id,
                     supplier_batch,
                     purity = -1):
    sSql = f'''insert into {bcpvsDB}.batch (
    compound_id,
    notebook_ref,
    suffix,
    submitter,
    submittal_date,
    project,
    purity,
    supplier,
    biological_mw,
    library_id,
    compound_type,
    product_type,
    supplier_id,
    supplier_batch,
    chemreg_regno)
    values (
    '{compound_id}',
    '{notebook_ref}',
    '{suffix}',
    '{submitter}',
    now(),
    '{project}',
    {purity},
    '{supplier}',
    {biological_mw},
    '{library_id}',
    '{compound_type}',
    '{product_type}',
    '{supplier_id}',
    '{supplier_batch}',
    {chemreg_regno})
    '''
    try:
        cur.execute(sSql)
    except Exception as e:
        logger.error(f"{sSql}")
        logger.error(f"{str(e)}")

    if suffix == '':
        sNULL = 'NULL'
        sSql = f"""update {bcpvsDB}.batch set
        suffix = {sNULL}
        where notebook_ref = '{notebook_ref}'"""
        cur.execute(sSql)


def getMoleculeProperties(self, molfile, chemregDB):
    sSql = f'''select
    bin2smiles(mol2bin('{molfile}'), 'mol') smiles,
    MolFormula(mol2bin('{molfile}', 'mol')),
    MolWeight(mol2bin(UNIQUEKEY('{molfile}', 'cistrans'))),
    MolNofMol(mol2bin('{molfile}', 'mol'))
    '''
    cur.execute(sSql)
    res = cur.fetchall()
    try:
        sSmiles = (res[0][0]).decode()
        C_MF = res[0][1].decode().replace(" ", "")
        mainFragMolWeight = res[0][2]
        iNrOfFragments = res[0][3]
    except:
        return (False, False, False, False, False, f'No smiles for molecule')

    try:
        molmassFormula = Formula(C_MF.replace('-', ''))
        C_CHNS = getAtomicComposition(molmassFormula.composition())
    except Exception as e:
        return (False, False, False, False, False, f'{str(e)}')

    if sSmiles == '':
        return (False, False, False, False, False, 'Empty molfile')
    if iNrOfFragments > 1:
        saSmileFragments = sSmiles.split('.')
    else:
        saSmileFragments = [sSmiles]

    if len(saSmileFragments) > 1 and iNrOfFragments == len(saSmileFragments):
        iFragPosition = 0
        for fragment in saSmileFragments:
            sSql = f'''select MolWeight(mol2bin(UNIQUEKEY('{fragment}', 'cistrans')))'''
            cur.execute(sSql)
            res = cur.fetchall()
            if res[0][0] == mainFragMolWeight:
                if iFragPosition != 0:
                    saSmileFragments[0], saSmileFragments[iFragPosition] = saSmileFragments[iFragPosition], saSmileFragments[0]
                    break
            iFragPosition += 1

    saSmileFragments = sorted(saSmileFragments, key=len, reverse=True)
    mainMolSmile = sSmiles
    saltSmile = ''

    if len(saSmileFragments) > 1:
        mainMolSmile = saSmileFragments[0]
        saltSmile = '.'.join(saSmileFragments[1:])

    cur.execute(f"""select
    MolWeight(mol2bin('{sSmiles}', 'smiles')),
    MolWeight(mol2bin('{mainMolSmile}', 'smiles'))
    """)
    resMolcart = cur.fetchall()
    C_MW = Formula(C_MF).mass
    C_MONOISO = resMolcart[0][1]
    if C_MONOISO == None:
        C_MONOISO = C_MW
    return (C_MF, C_MW, C_MONOISO, C_CHNS, saltSmile, '')


class PingDB(tornado.web.RequestHandler):
    def get(self):
        sSql = "select * from hive.user_details where pkey = 0"
        ret = cur.ping(sSql)
        if ret == 'error':
            self.set_status(400)

    def head(self):
        sSql = "select * from hive.user_details where pkey = 0"
        ret = cur.ping(sSql)
        if ret == 'error':
            self.set_status(400)


@jwtauth
class GetMolkeyct(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sCmpId = self.get_argument("compound_id")
        sSql = f'''
        select CONVERT(molkeyct USING utf8) from {bcpvsDB}.JCMOL_MOLTABLE_ukey
        where compound_id = '{sCmpId}'
        '''
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            res = res[0][0]
        else:
            res = ''
        self.finish(res)


@jwtauth
class GetCompoundDuplicates(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sCmpId = self.get_argument("compound_id")
        sSql = f'''
        select compound_id from {bcpvsDB}.JCMOL_MOLTABLE_ukey
        where molkeyct = (select molkeyct from {bcpvsDB}.JCMOL_MOLTABLE_ukey
                         where compound_id = '{sCmpId}')
        '''
        cur.execute(sSql)
        res = res2json()
        #res = cur.fetchall()
        if len(res) > 0:
            pass
            #res = res2json()
        else:
            res = ''
        self.finish(res)


@jwtauth
class CreateSalt(tornado.web.RequestHandler):
    def put(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sSmiles = self.get_argument("smiles")
        iSaltNumber = getNewSaltNumber(chemregDB)
        mol = Chem.MolFromSmiles(sSmiles)
        mw = Descriptors.MolWt(mol)
        mf = rdMolDescriptors.CalcMolFormula(mol)
        suffix = f'X{iSaltNumber}'
        sio = StringIO()
        with Chem.SDWriter(sio) as w:
            w.write(mol)
        molfile = sio.getvalue()
        sSql = f"""insert into {chemregDB}.salts (suffix, smiles, mw, mf, molfile)
        values ('{suffix}', '{sSmiles}', {mw}, '{mf}', '{molfile}')
        """
        cur.execute(sSql)


@jwtauth
class GetCanonicSmiles(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sSmiles = self.get_argument("smiles")
        sLetters, saRemainderSmile = getSaltLetters([sSmiles], chemregDB)
        mol = Chem.MolFromSmiles(sSmiles)
        flattenSmile = Chem.rdmolfiles.MolToSmiles(mol)
        canonSmile = Chem.CanonSmiles(flattenSmile)
        sRes = json.dumps([f'{sLetters}', canonSmile])
        self.finish(sRes)


@jwtauth
class GetLastBatchFromEln(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sEln = self.get_argument("eln")
        eln_id = sEln[0:6]
        sSql = f'''select notebook_ref from {bcpvsDB}.batch
        where notebook_ref like '{eln_id}%' order by notebook_ref desc'''
        cur.execute(sSql)
        cRes = cur.fetchall()
        if len(cRes) > 0:
            iBatch = int((cRes[0][0])[-3:])
            self.finish(f'{iBatch}')
        else:
            self.finish(b'0')


@jwtauth
class GetNextSdfSequence(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sSql = f"select pkey from {chemregDB}.sdfile_sequence"
        cur.execute(sSql)
        pkey = cur.fetchall()[0][0] +1
        sSql = f"update {chemregDB}.sdfile_sequence set pkey={pkey}"
        cur.execute(sSql)
        self.finish(f'{pkey}')


@jwtauth
class UpdateStructureAdmin(tornado.web.RequestHandler):
    def post(self):
        chemregDB, bcpvsDB = getDatabase(self)
        molfile = self.get_body_argument('molfile')
        compound_id = self.get_body_argument('compound_id')
        smiles = self.get_body_argument('smiles')
        database = bcpvsDB + '.JCMOL_MOLTABLE'
        (C_MF,
         C_MW,
         C_MONOISO,
         C_CHNS,
         saSalts,
         errorMessage) = getMoleculeProperties(self, molfile, chemregDB)
        sSql = f'''
        update {bcpvsDB}.compound
        set
        mf = '{C_MF}',
        sep_mol_monoiso_mass = '{C_MONOISO}',
        smiles_std = '{smiles}'
        where compound_id = '{compound_id}'
        '''
        cur.execute(sSql)

        sSql = f'''
update {database}
set mol = '{molfile}'
where compound_id = '{compound_id}'
        '''
        cur.execute(sSql)
        ######################################################
        ###  Update the molcart tables for compound

        sSql = f"""
update {database}_MOL
set mol = mol2bin('{molfile}', 'mol')
where compound_id = '{compound_id}'
        """
        cur.execute(sSql)

        sSql = f"""
update {database}_ukey
set
molkey = (select uniquekey(mol)
from {database} where compound_id = '{compound_id}'),
molkeyns = (select uniquekey(`mol`,'nostereo')
from {database} where compound_id = '{compound_id}'),
molkeyct = (select uniquekey(`mol`,'cistrans')
from {database} where compound_id = '{compound_id}')
where compound_id = '{compound_id}'
        """
        cur.execute(sSql)

        sSql = f"""
update {database}_MOL_keysim set molkey = (select fp(mol, 'sim')
from {database}_MOL where compound_id = '{compound_id}')
where compound_id = '{compound_id}'
        """
        cur.execute(sSql)

        sSql = f"""
update {database}_MOL_key set molkey = (select fp(mol, 'sss') as molkey
from {database}_MOL where compound_id = '{compound_id}')
where compound_id = '{compound_id}'
        """
        cur.execute(sSql)



@jwtauth
class GetRegnosFromSdfSequence(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sdfile_sequence = self.get_argument("sdfile_sequence")
        sSql = f'''select regno from {chemregDB}.chem_info
        where sdfile_sequence = {sdfile_sequence}
        '''
        cur.execute(sSql)
        res = res2json()
        self.finish(res)


@jwtauth
class BcpvsRegCompound(tornado.web.RequestHandler):
    def put(self):
        chemregDB, bcpvsDB = getDatabase(self)
        regno = self.get_argument("regno")
        sSql = f'''select regno,
        compound_id,
        c_mf,
        ip_rights,
        c_monoiso,
        molfile,
        jpage,
        suffix,
        chemist,
        project,
        source,
        c_mw,
        library_id,
        compound_type,
        product,
        external_id,
        supplier_batch,
        purity
        from {chemregDB}.chem_info
        where regno = {regno}'''
        cur.execute(sSql)
        (regno,
         compound_id,
         c_mf,
         ip_rights,
         c_monoiso,
         molfile,
         jpage,
         suffix,
         chemist,
         project,
         source,         #  Supplier
         c_mw,
         library_id,
         compound_type,
         product,
         external_id,
         supplier_batch,
         purity
        ) = cur.fetchall()[0]

        (C_MF,
         C_MW,
         C_MONOISO,
         C_CHNS,
         saSalts,
         errorMessage) = getMoleculeProperties(self, molfile, chemregDB)
        mol = Chem.MolFromMolBlock(molfile, removeHs=False, sanitize=True)

        molfile, stdSMILES = cleanStructureRDKit(molfile)
        if compound_id in ('', None):

            mols = checkUniqueStructure(stdSMILES, bcpvsDB)
            if mols == False:
                logger.error(f'Error in molfile for regno {self.regno}')
                return False

            if len(mols) != 0:
                compound_id = mols[0][0]
                compound_id_numeric = compound_id[3:]
            else:
                # Register a new compound
                compound_id_numeric = getNewCompoundId(bcpvsDB)

                compound_id = registerNewCompound(bcpvsDB,
                                                  compound_id_numeric,
                                                  molfile,
                                                  stdSMILES,
                                                  C_MF,
                                                  C_MONOISO,
                                                  ip_rights,
                                                  compound_name = '')
                # Is the moldepict on the next line the solution to stereo problems?
                #sSql = f'''select bin2mol(moldepict(mol2bin(UNIQUEKEY('{molfile}',
                #                                                      'cistrans'),
                #                                                      'smiles')))'''
                #sSql = f'''select bin2mol(mol2bin(UNIQUEKEY('{molfile}', 'cistrans'), 'smiles'))'''
                #cur.execute(sSql)
                #strippedMolfile = cur.fetchall()[0][0].decode("utf-8")
                strippedMolfile = molfile
                addStructure(f"{bcpvsDB}.JCMOL_MOLTABLE",
                             strippedMolfile,
                             compound_id,
                             'compound_id')
        sSql = f'''update {chemregDB}.chem_info
                   set compound_id='{compound_id}'
                   where regno = {regno}'''
        cur.execute(sSql)
        registerNewBatch(bcpvsDB,
                         compound_id,
                         regno,
                         jpage,
                         saSalts,
                         chemist,
                         project,
                         source,         #  Supplier
                         c_mw,
                         library_id,
                         compound_type,
                         product,
                         external_id,
                         supplier_batch,
                         purity = -1)
        self.finish(b'compound_id')


@jwtauth
class ChemRegAddMol(tornado.web.RequestHandler):
    def post(self):
        chemregDB, bcpvsDB = getDatabase(self)
        molfile = self.get_body_argument('molfile')
        jpage = self.get_body_argument('jpage')
        chemist = self.get_body_argument('chemist')
        compound_type = self.get_body_argument('compound_type')
        project = self.get_body_argument('project')
        source = self.get_body_argument('source')
        solvent = self.get_body_argument('solvent')
        product = self.get_body_argument('product')
        library_id = self.get_body_argument('library_id')
        sLib = re.search(r"(Lib-\d\d\d\d)", library_id)
        library_id = sLib.group()
        if library_id.startswith("Lib"):
            pass
        else:
            # If the library is blank set it to Compound collection lib
            library_id = 'Lib-3002'
        external_id = self.get_body_argument('external_id')
        supplier_batch = self.get_body_argument('supplier_batch')
        purity = self.get_body_argument('purity')
        ip_rights = self.get_body_argument('ip_rights')
        sdfile_sequence = self.get_body_argument('sdfile_sequence')

        if "0  0  0     0  0            999 V2" in molfile:
            self.set_status(500)
            self.finish(f'Nostruct for jpage: {jpage}')
            return

        try:
            l1, l2= cleanStructureRDKit(molfile)
            if l1 == False:
                self.set_status(500)
                self.finish(f'RDKit failed to convert Molfile for: {jpage}')
                return                
        except:
            self.set_status(500)
            self.finish(f'RDKit failed to convert Molfile for: {jpage}')
            return
            
            
        (C_MF,
         C_MW,
         C_MONOISO,
         C_CHNS,
         saSalts,
         errorMessage) = getMoleculeProperties(self, molfile, chemregDB)
        if saSalts == None:
            saSalts = ''
        if C_MF == False:
            self.set_status(500)
            self.finish(f'Molfile failed {external_id} {errorMessage}')
            logger.error(f'Molfile failed {external_id} {errorMessage}')
            return
        if saSalts == False:
            self.set_status(500)
            self.finish(f'Unknown salt in molfile for {external_id} {errorMessage}')
            logger.error(f'Unknown salt in molfile for {external_id} {errorMessage}')
            return

        ########################
        # Do exact match with molecule against the CBK database

        #sSql = f'''select bin2mol(moldepict(mol2bin(UNIQUEKEY('{molfile}', 'cistrans'), 'smiles')))'''
        #cur.execute(sSql)
        #strip = cur.fetchall()[0][0].decode("utf-8")

        molfileRdkit, stdSMILES = cleanStructureRDKit(molfile)
        mols = checkUniqueStructure(stdSMILES, bcpvsDB)
        #mols = checkUniqueStructure(strip, bcpvsDB)
        #mols = checkUniqueStructure(molfile, bcpvsDB)
        if mols == False:
            return False
        sStatus = 'oldMolecule'
        if len(mols) == 0:
            sStatus = 'newMolecule'
            compound_id = ''
        else:
            compound_id = mols[0][0]

        ########################
        # Get new regno and add nmolecule to chem_info
        newRegno = getNewRegno(chemregDB)
        if purity == '':
            purity = -1
        sSql = f"""
        insert into {chemregDB}.chem_info (
        regno,
        jpage,
        compound_id,
        rdate,
        chemist,
        compound_type,
        project,
        source,
        solvent,
        product,
        library_id,
        external_id,
        supplier_batch,
        purity,
        ip_rights,
        C_CHNS,
        C_MF,
        C_MW,
        C_MONOISO,
        SUFFIX,
        sdfile_sequence,
        molfile)
        values (
        '{newRegno}',
        '{jpage}',
        '{compound_id}',
        now(),
        '{chemist}',
        '{compound_type}',
        '{project}',
        '{source}',
        '{solvent}',
        '{product}',
        '{library_id}',
        '{external_id}',
        '{supplier_batch}',
        {purity},
        '{ip_rights}',
        '{C_CHNS}',
        '{C_MF}',
        {C_MW},
        {C_MONOISO},
        '{saSalts}',
        {sdfile_sequence},
        '{molfile}'
        )
        """
        cur.execute(sSql)
        if saSalts == '':
            sNULL = 'NULL'
            sSql = f"""update {chemregDB}.chem_info set
            suffix = {sNULL}
            where regno = {newRegno}"""
            cur.execute(sSql)

        addStructure(f"{chemregDB}.CHEM", molfile, newRegno, 'regno')
        self.finish(sStatus)


@jwtauth
class Search(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        column = self.get_argument("column")
        value = self.get_argument("value")
        values = (value, )

        if column == 'JPAGE':
            sSql = f"select regno from {chemregDB}.chem_info where {column} like '{value}%'"
        else:
            sSql = f"select regno from {chemregDB}.chem_info where {column} = '{value}'"
        cur.execute(sSql)
        res = res2json()
        self.write(res)


@jwtauth
class GetRegnoData(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        column = self.get_argument("column")
        regno = self.get_argument("regno")
        values = (regno, )
        sSql = f"select {column} from {chemregDB}.chem_info where regno = {regno}"
        cur.execute(sSql)
        res = res2json()
        self.write(res)


@jwtauth
class CreateMolImageFromMolfile(tornado.web.RequestHandler):
    def post(self):
        sStatus = ''
        molfile = '\n' + self.get_body_argument('molfile')
        sMolId = self.get_body_argument('mol_id')
        try:
            os.remove(f'mols/{sMolId}.png')
        except:
            pass
        if len(molfile) < 85:
            # The molecule is a nostruct, rdkit doesn't like that
            resDict = {
                "molfile": '',
                "status": '',
                "smiles": ''
            }
            self.finish(json.dumps(resDict))
            return
        mol = Chem.MolFromMolBlock(molfile, removeHs=False, sanitize=True)
        params = rdMolStandardize.CleanupParameters()
        params.tautomerRemoveSp3Stereo = False
        params.tautomerRemoveBondStereo = False
        params.tautomerRemoveIsotopicHs = False
        try:
            clean_mol = rdMolStandardize.Cleanup(mol, params)
        except:
            logger.error(f'Failed to standardize molfile {sMolId}')
            sStatus = 'Molecule altered, checkit'
            sSql = f'''select bin2mol(moldepict(mol2bin(UNIQUEKEY('{molfile}',
                                                                  'cistrans'),
                                                                  'smiles')))'''
            cur.execute(sSql)
            molfile = cur.fetchall()[0][0].decode("utf-8")
            mol = Chem.MolFromMolBlock(molfile, removeHs=False, sanitize=True)
            if mol == None:
                logger.error(f'Failed to standardize molfile: {molfile}')
                return
            params = rdMolStandardize.CleanupParameters()
            params.tautomerRemoveSp3Stereo = False
            params.tautomerRemoveBondStereo = False
            params.tautomerRemoveIsotopicHs = False
            clean_mol = rdMolStandardize.Cleanup(mol, params)

        parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol, params)
        uncharger = rdMolStandardize.Uncharger()
        uncharged_parent_clean_mol = uncharger.uncharge(parent_clean_mol)
        te = rdMolStandardize.TautomerEnumerator(params) # idem
        taut_uncharged_parent_clean_mol = te.Canonicalize(uncharged_parent_clean_mol)
        smiles = Chem.MolToSmiles(taut_uncharged_parent_clean_mol)
        m = Chem.MolToMolBlock(taut_uncharged_parent_clean_mol)

        # Don't think we need to generate new coords with molcart
        #sSql = f'''select bin2mol(moldepict(mol2bin(UNIQUEKEY('{m}', 'cistrans'), 'smiles')))'''
        #cur.execute(sSql)
        #strip = cur.fetchall()[0][0].decode("utf-8")
        #m = Chem.MolFromMolBlock(strip)

        m = 'Id' + m
        strip = m
        new_mol = Chem.MolFromMolBlock(strip, removeHs=False, sanitize=True)
        try:
            Draw.MolToFile(new_mol, f'mols/{sMolId}.png', kekulize=True, size=(280, 280))
        except Exception as e:
            logger.error(str(e))
            logger.error(f"{sMolId} is nostruct, png failed")
        resDict = {
            "molfile": strip,
            "status": sStatus,
            "smiles": smiles
        }

        self.finish(json.dumps(resDict))

@jwtauth
class CreateMolImage(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        regno = self.get_argument("regno")
        sSql = f"""select molfile from {chemregDB}.chem_info
                   where regno = '{regno}'"""
        cur.execute(sSql)
        molfile = cur.fetchall()

        if len(molfile) == 0:
            sSql = f"""select molfile from chemspec.chem_info
            where regno = '{regno}'"""
            cur.execute(sSql)
            molfile = cur.fetchall()

        if len(molfile) > 0 and molfile[0][0] != None:
            try:
                createPngFromMolfile(regno, molfile[0][0])
            except:
                pass
        self.finish()


@jwtauth
class CreateBcpvsMolImage(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sCmpId = self.get_argument("compound_id")
        sSql = f"""select mol, compound_id from {bcpvsDB}.JCMOL_MOLTABLE
                   where compound_id like '{sCmpId}%'"""
        cur.execute(sSql)
        molfile = cur.fetchall()
        if len(molfile) > 0 and molfile[0][0] != None:
            try:
                createPngFromMolfile(sCmpId, molfile[0][0])
            except:
                pass
        else:
            try:
                os.remove(f'mols/{sCmpId}.png')
            except:
                pass
        self.finish()


@jwtauth
class LoadMolfile(tornado.web.RequestHandler):
    def post(self):
        chemregDB, bcpvsDB = getDatabase(self)
        fBody = self.request.files['file'][0]
        regnoBody = self.request.files['regno'][0]
        molfile = tornado.escape.xhtml_unescape(fBody.body)
        regno = tornado.escape.xhtml_unescape(regnoBody.body)

        (C_MF,
         C_MW,
         C_MONOISO,
         C_CHNS,
         saSalts,
         errorMessage) = getMoleculeProperties(self, molfile, chemregDB)
        if C_MF == False:
            self.set_status(500)
            self.finish(f'Molfile failed {regno} {errorMessage}')
            return
        try:
            addStructure(f'{chemregDB}.CHEM', molfile, regno, 'regno')
        except Exception as e:
            logger.error(str(e))
            logger.error(f'failed on regno {regno}')
            return
        compound_id = isItNewStructure(self, molfile)
        sSql = f"""update {chemregDB}.chem_info set
                   `molfile` = '{molfile}',
                    compound_id = '{compound_id}',
                    C_MF = '{C_MF}',
                    C_MW = {C_MW},
                    C_MONOISO = {C_MONOISO},
                    C_CHNS = '{C_CHNS}',
                    suffix = '{saSalts}'
                  where regno = {regno}"""
        cur.execute(sSql)

        if saSalts == '':
            sNULL = 'NULL'
            sSql = f"""update {chemregDB}.chem_info set
            suffix = {sNULL}
            where regno = {regno}"""
            cur.execute(sSql)

        createPngFromMolfile(regno, molfile)


@jwtauth
class UpdateRegnoBatch(tornado.web.RequestHandler):
    def put(self):
        chemregDB, bcpvsDB = getDatabase(self)
        regno = self.get_argument("regno")
        batch = self.get_argument("batch")
        sSql = f"""select compound_id from {bcpvsDB}.batch
                   where notebook_ref = '{batch}'"""
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            self.set_status(500)
            self.finish()
        else:
            sSql = f"""select regno from {chemregDB}.chem_info
            where jpage = '{batch}' and regno != '{regno}'"""
            cur.execute(sSql)
            res = cur.fetchall()
            if len(res) > 0:
                self.set_status(500)
                self.finish()
            else:
                sSql = f"""update {chemregDB}.chem_info
                set jpage = '{batch}' where regno = '{regno}'"""
                cur.execute(sSql)
                self.finish()


@jwtauth
class CreateSupplier(tornado.web.RequestHandler):
    def put(self):
        chemregDB, bcpvsDB = getDatabase(self)
        supplierName = self.get_argument("supplier")
        sSql = f"""insert into {bcpvsDB}.compound_suppliers
        (name) values ('{supplierName}')"""
        cur.execute(sSql)


@jwtauth
class CreateRegno(tornado.web.RequestHandler):
    def put(self):
        chemregDB, bcpvsDB = getDatabase(self)
        regno = self.get_argument("regno")
        sSql = f"""insert into {chemregDB}.chem_info
        (regno, rdate) values (%s, now())"""
        val = (regno, )
        cur.execute(sSql, val)


@jwtauth
class DeleteRegno(tornado.web.RequestHandler):
    def put(self):
        chemregDB, bcpvsDB = getDatabase(self)
        regno = self.get_argument("regno")
        sSql = f"""select regno, c_mf, chemist from {chemregDB}.chem_info
                  where c_mf is null and chemist is null and
                  regno={regno}"""
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            logger.info('Deleting ' + str(res))
            sSql = f"""delete from {chemregDB}.chem_info
                      where regno = {regno}"""
            cur.execute(sSql)


@jwtauth
class UpdateColumn(tornado.web.RequestHandler):
    def put(self):
        chemregDB, bcpvsDB = getDatabase(self)
        column = self.get_argument("column")
        value = self.get_argument("value")
        regno = self.get_argument("regno")

        if 'PURITY' in column:
            try:
                value = int(value)
            except:
                logger.info(f"Purity is non integer, setting purity to -1")
                value = -1
        values = (value, regno, )
        sSql = f"""update {chemregDB}.chem_info set {column} = '{value}'
        where regno = {regno}"""

        if column == 'library_id':
            sLib = re.search(r"(Lib-\d\d\d\d)", value)
            try:
                value = sLib.group()
                sSql = f"""update {chemregDB}.chem_info set {column} = '{value}'
                where regno = {regno}"""
            except:
                sSql = f"""update {chemregDB}.chem_info set {column} = NULL
                where regno = {regno}"""

        cur.execute(sSql)


@jwtauth
class RegCompound(tornado.web.RequestHandler):
    def get(self):
        # Contains user found in previous auth
        if self.request.headers.get('auth'):
            self.write('ok')


@jwtauth
class GetCompound(tornado.web.RequestHandler):
    def get(self):
        # Contains user found in previous auth
        if self.request.headers.get('auth'):
            self.write('ok')


@jwtauth
class GetMolfile(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        regno = self.get_argument("regno")
        sSql = f"""select molfile from {chemregDB}.chem_info
                   where regno = '{regno}'"""
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            self.write(res[0][0])
            return

        sSql = f"""select molfile from chemspec.chem_info
                   where regno = '{regno}'"""
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            self.write(res[0][0])


@jwtauth
class GetForwardCompound(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sCmpId = self.get_argument("compound_id")
        sSql = f"""select compound_id from {bcpvsDB}.compound
                   where compound_id > '{sCmpId}' order by compound_id asc limit 1"""
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            self.write(json.dumps(res[0][0]))


@jwtauth
class GetBackwardCompound(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sCmpId = self.get_argument("compound_id")
        sSql = f"""select compound_id from {bcpvsDB}.compound
                   where compound_id < '{sCmpId}' order by compound_id desc limit 1"""
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            self.write(json.dumps(res[0][0]))


@jwtauth
class GetBackwardRegno(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        regno = self.get_argument("regno")

        sSql = f"""select greatest(chemreg, chemspec) from (
        WITH
        cte1 AS (SELECT IFNULL((SELECT regno FROM {chemregDB}.chem_info
        where regno < '{regno}' order by regno desc limit 1), 0) as chemreg),
        cte2 AS (select IFNULL((SELECT regno chemspec FROM chemspec.chem_info
        where regno < '{regno}' order by regno desc limit 1), 0) as chemspec)
        SELECT chemreg, chemspec FROM cte1, cte2) ss
        """
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            self.write(json.dumps(res[0][0]))


@jwtauth
class GetForwardRegno(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        regno = self.get_argument("regno")

        sSql = f"""select least(chemreg, chemspec) from (
        WITH
        cte1 AS (SELECT regno chemreg FROM {chemregDB}.chem_info
        where regno > '{regno}' order by regno asc limit 1),
        cte2 AS (select IFNULL((SELECT regno chemspec FROM chemspec.chem_info
        where regno > '{regno}' order by regno asc limit 1), 999999999) as chemspec)
        SELECT chemreg, chemspec FROM cte1, cte2) ss
        """
        cur.execute(sSql)
        res = cur.fetchall()

        if len(res) > 0:
            self.write(json.dumps(res[0][0]))


@jwtauth
class GetRegnoFromCompound(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sCmpId = self.get_argument("compound_id")
        if sCmpId.startswith('CBK6'):
            sSql = f"""select regno from {chemregDB}.chem_info
                       where compound_id = '{sCmpId}'"""
        else:
            sSql = f"""
            select cs.regno from chemspec.chem_info cs
            join {bcpvsDB}.batch b on cs.regno = b.chemspec_regno
            where b.compound_id = '{sCmpId}'
                    """

        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            self.write(json.dumps(res[0][0]))


@jwtauth
class GetCompoundFromRegno(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sRegno = self.get_argument("regno")
        sSql = f"""select compound_id from {chemregDB}.chem_info
                   where regno = '{sRegno}'"""
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            self.write(json.dumps(res[0][0]))
            return

        sSql = f"""
        select compound_id from {bcpvsDB}.batch b
        join chemspec.chem_info ci on ci.regno = b.chemspec_regno
        where regno = '{sRegno}'"""
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            self.write(json.dumps(res[0][0]))


@jwtauth
class GetMolfileBcpvs(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        compound_id = self.get_argument("compound_id")
        sSql = f"""select mol from {bcpvsDB}.JCMOL_MOLTABLE
                   where compound_id = '{compound_id}'"""
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            self.write(res[0][0])


@jwtauth
class GetTextColumn(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        column = self.get_argument("column")
        regno = self.get_argument("regno")
        if regno == None:
            return

        if column == 'library_id':
            sSql = f"""
select IFNULL((select concat(library_id, " ", l.description)
from {chemregDB}.chem_info c, {bcpvsDB}.compound_library l
where c.library_id = l.library_name
and c.regno = '{regno}'), null)
            """
        else:
            sSql = f"""select cast({column} as char)
            from {chemregDB}.chem_info where regno = {regno}"""
        cur.execute(sSql)
        res = res2json()
        self.write(res)


@jwtauth
class GetColComboData(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        column = self.get_argument("column")
        sSql = ''
        if column == 'project':
            sSql = """select project_name from hive.project_details
                      order by project_name"""
            sSql = """select project_name from hive.project_details
                      order by created_date desc"""
        elif column == 'chemist':
            sSql = """select fullname from hive.user_details
            where ORGANIZATION = 'chemistry'"""
        elif column == 'compound_type':
            sSql = f"select type from {bcpvsDB}.compound_type order by type"
        elif column == 'product':
            sSql = f"SELECT type FROM {bcpvsDB}.product_type order by type"
        elif column == 'supplier':
            sSql = f"SELECT name FROM {bcpvsDB}.compound_suppliers order by name"
        elif column == 'solvent':
            sSql = "SELECT solvent FROM chemspec.solvent_tbl order by solvent"
        elif column == 'library_id':
            sSql = f"""SELECT concat(library_name, " ", description) as library_name
                       FROM {bcpvsDB}.compound_library order by library_name"""
        elif column == 'library_description':
            sSql = f"""SELECT description FROM {bcpvsDB}.compound_library
                      order by description"""
        cur.execute(sSql)
        res = res2json()
        self.write(res)


@jwtauth
class GetLibraryName(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        try:
            library_id = self.get_argument("library_id")
        except:
            library_id = ''
        sLib = re.search(r"(Lib-\d\d\d\d)", library_id)
        if sLib == None:
            self.write(b'[]')
            return
        library_id = sLib.group()
        sSql = f"""SELECT description FROM {bcpvsDB}.compound_library
                  where library_name = '{library_id}'"""
        cur.execute(sSql)
        res = res2json()
        self.write(res)


@jwtauth
class GetLibraries(tornado.web.RequestHandler):
    def get(self):
        sSql = "select fullname from hive.user_details where ORGANIZATION = 'chemistry'"
        cur.execute(sSql)
        res = res2json()
        self.write(res)


@jwtauth
class GetNextRegno(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        newRegno = getNewRegno(chemregDB)
        self.write(json.dumps(newRegno))


@jwtauth
class CreateLibrary(tornado.web.RequestHandler):
    def put(self):
        chemregDB, bcpvsDB = getDatabase(self)
        sLibraryDescription = self.get_argument("library_name")
        sSupplier = self.get_argument("supplier")
        restricted = self.get_argument("restricted")
        if restricted == 'False':
            restricted = 'No'
        else:
            restricted = 'Yes'
            
        sSql = f"select pkey from {bcpvsDB}.library_id_sequence"
        cur.execute(sSql)
        pkey = cur.fetchall()[0][0] +1
        sSql = f"update {bcpvsDB}.library_id_sequence set pkey={pkey}"
        cur.execute(sSql)

        library_id = f'Lib-{pkey}'

        sSql = f"""insert into {bcpvsDB}.compound_library
        (library_name,
        supplier,
        description,
        restricted) values
        ('{library_id}',
        '{sSupplier}',
        '{sLibraryDescription}',
        '{restricted}')
        """
        cur.execute(sSql)


class GetDatabase(tornado.web.RequestHandler):
    def get(self):
        sRes = json.dumps([['Live'], ['Test']])
        self.finish(sRes)

@jwtauth
class UploadBinary(tornado.web.RequestHandler):
    def post(self, *args, **kwargs):
        os_name = self.get_argument("os_name")

        try:
            # self.request.files['file'][0]:
            # {'body': 'Label Automator ___', 'content_type': u'text/plain', 'filename': u'k.txt'}
            file1 = self.request.files['file'][0]
        except:
            logging.error("Error cant find file1 in the argument list")
            return

        bin_file = ""
        if os_name == 'Windows':
            bin_file = f'dist/{os_name}/ch.exe'
        elif os_name == 'Linux':
            bin_file = f'dist/{os_name}/ch'
        elif os_name == 'Darwin':
            bin_file = f'dist/{os_name}/ch'
        else:
            # unsupported OS
            self.set_status(500)
            self.write({'message': 'OS not supported'})
            return

        output_file = open(bin_file, 'wb')
        output_file.write(file1['body'])
        output_file.close()

class GetChemRegBinary(tornado.web.RequestHandler):
    def post(self):
        pass

    def get(self, os_name):
        bin_file = ""
        if os_name == 'Windows':
            bin_file = f'dist/{os_name}/ch.exe'
        elif os_name == 'Linux':
            bin_file = f'dist/{os_name}/ch'
        elif os_name == 'Darwin':
            bin_file = f'dist/{os_name}/ch'
        else:
            # unsupported OS
            self.set_status(500)
            self.write({'message': 'OS not supported'})
            return
        try:
            with open(bin_file, 'rb') as f:
                logging.info("sending bin file")
                self.set_status(200)
                self.write(f.read())
        except Exception as e:
            logging.error(f"Did not send bin file, error: {str(e)}")

class getVersionData(tornado.web.RequestHandler):
    def post(self):
        pass

    def get(self):
        try:
            with open('./ver.dat', 'r') as f:
                self.write(json.load(f))
                return
        except Exception as e:
            logging.error(str(e))
            self.set_status(500)
            self.write({'message': 'ver.dat not available'})

@jwtauth
class UploadVersionNo(tornado.web.RequestHandler):
    def post(self, *args, **kwargs):
        ver_no = self.get_argument("ver_no")
        ver_file = "ver.dat"

        with open(ver_file, "r") as f:
            data = json.load(f)

        data["version"] = ver_no

        with open(ver_file, "w") as f:
            json.dump(data, f)

@jwtauth
class UploadLauncher(tornado.web.RequestHandler):
    def post(self, *args, **kwargs):
        os_name = self.get_argument("os_name")

        try:
            # self.request.files['file'][0]:
            # {'body': 'Label Automator ___', 'content_type': u'text/plain', 'filename': u'k.txt'}
            file1 = self.request.files['file'][0]
        except:
            logging.error("Error cant find file1 in the argument list")
            return

        bin_file = ""
        if os_name == 'Windows':
            bin_file = f'dist/launchers/{os_name}/chemreg.exe'
        elif os_name == 'Linux':
            bin_file = f'dist/launchers/{os_name}/chemreg'
        elif os_name == 'Darwin':
            bin_file = f'dist/launchers/{os_name}/chemreg'
        else:
            # unsupported OS
            self.set_status(500)
            self.write({'message': 'OS not supported'})
            return

        output_file = open(bin_file, 'wb')
        output_file.write(file1['body'])
        output_file.close()


class ChemblExport(tornado.web.RequestHandler):
    def get(self, sRIDX, sBatches):
        pass
    def post(self, *args, **kwargs):
        try:
            sRIDX = self.get_argument("RIDX").strip()
            sAIDX = self.get_argument("AIDX").strip()
            sProject = self.get_argument("project").strip()
            sELN = self.get_argument("ELN").strip()
            sTARGET_TYPE = self.get_argument("TARGET_TYPE").strip()
            sASSAY_TYPE = self.get_argument("ASSAY_TYPE").strip()
            sACTION_TYPE = self.get_argument("ACTION_TYPE").strip()
            sBatches = self.get_argument("batches").strip()
        except:
            logging.error(f"chembl error")
            return
        '''
        if len(sProject) > 2 and len(sELN) > 2:
            sBatches = ''
        '''
        saBatches = sBatches.split()
        sMol = ''
        random_number = str(random.randint(0, 1000000))
        dir_name = f'dist/export/{random_number}'
        os.makedirs(dir_name, exist_ok=True)
        file_path = dir_name + "/COMPOUND_RECORD.tsv"
        molfile_path = dir_name + "/COMPOUND_CTAB.sdf"
        
        with open(file_path, 'w') as compound_record_file, open(molfile_path, 'w') as molfile_file:
            sHeader = ('CIDX', 'RIDX', 'COMPOUND_NAME', 'COMPOUND_KEY')

            compound_record_file.write('\t'.join(map(str, sHeader)) + '\n')
            if len(saBatches) > 10000000:
                chembl_export.exportFromBatches(cur, saBatches, sRIDX, compound_record_file, molfile_file)
            elif len(sProject) > 2 and len(sELN) > 2:
                saCompounds = saBatches
                chembl_export.exportFromElnProject(cur,
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
                                                   logging)

        assay_tsv_name = dir_name + "/ASSAY.tsv"
        with open(assay_tsv_name, 'w') as assay_tsv_file:
            sHeader = (
                'AIDX',
                'RIDX',
                'ASSAY_DESCRIPTION',
                'ASSAY_TYPE',
                'ASSAY_ORGANISM',
                'ASSAY_STRAIN',
                'ASSAY_TAX_ID',
                'ASSAY_TISSUE',
                'ASSAY_CELL_TYPE',
                'ASSAY_SUBCELLULAR_FRACTION',
                'ASSAY_SOURCE',
                'TARGET_TYPE',
                'TARGET_NAME',
                'TARGET_ACCESSION',
                'TARGET_ORGANISM',
                'TARGET_TAX_ID')
            # sAIDX =
            # sRIDX =
            sASSAY_DESCRIPTION = ''
            # sASSAY_TYPE =
            sASSAY_ORGANISM = ''
            sASSAY_STRAIN = ''
            sASSAY_TAX_ID = ''
            sASSAY_TISSUE = ''
            sASSAY_CELL_TYPE = ''
            sASSAY_SUBCELLULAR_FRACTION = ''
            sASSAY_SOURCE = ''
            # sTARGET_TYPE = 
            sTARGET_NAME = ''
            sTARGET_ACCESSION = ''
            sTARGET_ORGANISM = ''
            sTARGET_TAX_ID = ''
            assay_tsv_file.write('\t'.join(map(str, sHeader)) + '\n')
            assay_tsv_file.write(f'''{sAIDX}\t{sRIDX}\t{sASSAY_DESCRIPTION}\t{sASSAY_TYPE}\t{sASSAY_ORGANISM}\t{sASSAY_STRAIN}\t{sASSAY_TAX_ID}\t{sASSAY_TISSUE}\t{sASSAY_CELL_TYPE}\t{sASSAY_SUBCELLULAR_FRACTION}\t{sASSAY_SOURCE}\t{sTARGET_TYPE}\t{sTARGET_NAME}\t{sTARGET_ACCESSION}\t{sTARGET_ORGANISM}\t{sTARGET_TAX_ID}\n''')
        

        reference_tsv_name = dir_name + "/REFERENCE.tsv"
        with open(reference_tsv_name, 'w') as reference_tsv_file:
            sHeader = (
                'RIDX',
                'PUBMED_ID',
                'JOURNAL_NAME',
                'YEAR',
                'VOLUME',
                'ISSUE',
                'FIRST_PAGE',
                'LAST_PAGE',
                'REF_TYPE',
                'TITLE',
                'DOI',
                'PATENT_ID',
                'ABSTRACT',
                'AUTHORS')
            reference_tsv_file.write('\t'.join(map(str, sHeader)) + '\n')
            reference_tsv_file.write(f'''{sRIDX}\t\t\t\t\t\t\t\t\t\t\t\t\t\n''')


        assay_param_tsv_name = dir_name + "/ASSAY_PARAM.tsv"
        with open(assay_param_tsv_name, 'w') as assay_param_tsv_file:
            sHeader = (
                'AIDX',
                'TYPE',
                'RELATION',
                'VALUE',
                'UNITS',
                'TEXT_VALUE',
                'COMMENTS')

            assay_param_tsv_file.write('\t'.join(map(str, sHeader)) + '\n')
            assay_param_tsv_file.write(f'''{sAIDX}\t\t\t\t\t\t\n''')

            
        zip_filepath = os.path.join(dir_name, "ChEMBL_deposit.zip")

        with zipfile.ZipFile(zip_filepath, 'w') as zipf:
            zipf.write(os.path.join(dir_name, "COMPOUND_RECORD.tsv"), arcname="COMPOUND_RECORD.tsv")
            zipf.write(os.path.join(dir_name, "COMPOUND_CTAB.sdf"), arcname="COMPOUND_CTAB.sdf")
            zipf.write(os.path.join(dir_name, "ACTIVITY.tsv"), arcname="ACTIVITY.tsv")
            zipf.write(os.path.join(dir_name, "ASSAY.tsv"), arcname="ASSAY.tsv")
            zipf.write(os.path.join(dir_name, "REFERENCE.tsv"), arcname="REFERENCE.tsv")
            zipf.write(os.path.join(dir_name, "ASSAY_PARAM.tsv"), arcname="ASSAY_PARAM.tsv")

        
        self.write(json.dumps(f'''<a href=https://esox3.scilifelab.se/chemreg/dist/export/{random_number}/ChEMBL_deposit.zip>ChEMBL_deposit.zip</a>'''))
