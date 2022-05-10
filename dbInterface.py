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
from molmass import Formula
from io import StringIO
import sys
import codecs
import re

logger = logging.getLogger(__name__)

db_connection = MySQLdb.connect(
    host=config.database['host'],
    user=config.database['user'],
    passwd=config.database['password']
)
db_connection.autocommit(True)
cur = db_connection.cursor()

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

def checkUniqueStructure(molfile, bcpvsDB):
    sSql = f"""
SELECT * FROM {bcpvsDB}.`JCMOL_MOLTABLE_ukey` T1
WHERE T1.molkeyct = UNIQUEKEY('{molfile}', 'cistrans')
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
    m = Chem.MolFromMolBlock(molfile)
    try:
        Draw.MolToFile(m, f'mols/{regno}.png', size=(280, 280))
    except:
        logger.error(f"regno {regno} is nostruct")

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
    sep_mol_monoiso_mass)
    values (
    '{compound_id}',
    {compound_id_numeric},
    now(),
    '{mf}',
    '{ip_rights}',
    {sep_mol_monoiso_mass})
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
    chemspec_regno)
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
            sSql = f'''select MolWeight(mol2bin(UNIQUEKEY('{fragment}', 'cistrans')))
            '''
            cur.execute(sSql)
            res = cur.fetchall()
            if res[0][0] == mainFragMolWeight:
                if iFragPosition != 0:
                    saSmileFragments[0], saSmileFragments[iFragPosition] = saSmileFragments[iFragPosition], saSmileFragments[0]
                    break
            iFragPosition += 1

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
    C_MW = resMolcart[0][0]
    C_MONOISO = resMolcart[0][1]
    if C_MONOISO == None:
        C_MONOISO = C_MW
    return (C_MF, C_MW, C_MONOISO, C_CHNS, saltSmile, '')


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
        mol = Chem.MolFromMolBlock(molfile)
        if compound_id in ('', None):
            mols = checkUniqueStructure(molfile, bcpvsDB)
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
                                                  C_MF,
                                                  C_MONOISO,
                                                  ip_rights,
                                                  compound_name = '')
                addStructure(f"{bcpvsDB}.JCMOL_MOLTABLE",
                             molfile,
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
        external_id = self.get_body_argument('external_id')
        supplier_batch = self.get_body_argument('supplier_batch')
        purity = self.get_body_argument('purity')
        ip_rights = self.get_body_argument('ip_rights')
        sdfile_sequence = self.get_body_argument('sdfile_sequence')

        if "0  0  0     0  0            999 V2" in molfile:
            self.set_status(500)
            self.finish(f'Nostruct for jpage: {jpage}')
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
        mols = checkUniqueStructure(molfile, bcpvsDB)
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
        addStructure(f"{chemregDB}.CHEM", molfile, newRegno, 'regno')
        self.finish(sStatus)
          

@jwtauth
class Search(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        column = self.get_argument("column")
        value = self.get_argument("value")
        values = (value, )

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
class CreateMolImage(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        regno = self.get_argument("regno")
        sSql = f"""select molfile from {chemregDB}.chem_info
                   where regno = '{regno}'"""
        cur.execute(sSql)
        molfile = cur.fetchall()
        if len(molfile) > 0 and molfile[0][0] != None:
            createPngFromMolfile(regno, molfile[0][0])
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
        sSql = f"""update {chemregDB}.chem_info set
                   `molfile` = '{molfile}',
                    C_MF = '{C_MF}',
                    C_MW = {C_MW},
                    C_MONOISO = {C_MONOISO},
                    C_CHNS = '{C_CHNS}',
                    suffix = '{saSalts}'
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
        values = (value, regno, )
        if column == 'library_id':
            sLib = re.search(r"(Lib-\d\d\d\d)", value)
            value = sLib.group()

        sSql = f"""update {chemregDB}.chem_info set {column} = '{value}'
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


@jwtauth
class GetTextColumn(tornado.web.RequestHandler):
    def get(self):
        chemregDB, bcpvsDB = getDatabase(self)
        column = self.get_argument("column")
        regno = self.get_argument("regno")

        if column == 'library_id':
            sSql = f"""
select concat(library_id, " ", l.description)
from {chemregDB}.chem_info c, {bcpvsDB}.compound_library l
where c.library_id = l.library_name
and c.regno = '{regno}'
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
        library_id = self.get_argument("library_id")
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

        sSql = f"select pkey from {bcpvsDB}.library_id_sequence"
        cur.execute(sSql)
        pkey = cur.fetchall()[0][0] +1
        sSql = f"update {bcpvsDB}.library_id_sequence set pkey={pkey}"
        cur.execute(sSql)

        library_id = f'Lib-{pkey}'
        
        sSql = f"""insert into {bcpvsDB}.compound_library
        (library_name,
        supplier,
        description) values
        ('{library_id}',
        '{sSupplier}',
        '{sLibraryDescription}')
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
