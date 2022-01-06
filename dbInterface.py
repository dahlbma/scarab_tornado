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

logger = logging.getLogger(__name__)

db_connection = MySQLdb.connect(
    host=config.database['host'],
    user=config.database['user'],
    passwd=config.database['password'],
    database='chem_reg'
)
db_connection.autocommit(True)
cur = db_connection.cursor()

#sSql = f'select pkey, suffix, smiles, mf, mw from chem_reg.salts order by pkey'
#cur.execute(sSql)
#salts = cur.fetchall()

def sqlExec(sSql, values=None):
    if values == None:
        cur.execute(sSql)
    else:
        cur.execute(sSql, values)
    result = [list(i) for i in cur.fetchall()]
    return json.dumps(result)

def checkUniqueStructure(smiles):
    sSql = f"""
    select compound_id, bin2smiles(bcpvs.jcmol_moltable.mol) from
    bcpvs.jcmol_moltable_ukey join bcpvs.jcmol_moltable on
    (bcpvs.jcmol_moltable.molid=bcpvs.jcmol_moltable_ukey.molid)
    where uniquekey(mol2bin('{smiles}', 'smiles'))=molkey
    """
    cur.execute(sSql)
    mols = cur.fetchall()
    return mols

def getNewRegno():
    sSql = "select pkey from chem_reg.regno_sequence"
    cur.execute(sSql)
    pkey = cur.fetchall()[0][0] +1
    sSql = f"update chem_reg.regno_sequence set pkey={pkey}"
    cur.execute(sSql)
    return pkey

def getNewSaltNumber():
    sSql = "select pkey from chem_reg.salts_sequence"
    cur.execute(sSql)
    pkey = cur.fetchall()[0][0] +1
    sSql = f"update chem_reg.salts_sequence set pkey={pkey}"
    cur.execute(sSql)
    return pkey

def getNewCompoundId():
    sSql = "select pkey from bcpvs.compound_id_sequence"
    cur.execute(sSql)
    pkey = cur.fetchall()[0][0] +1
    sSql = f"update bcpvs.compound_id_sequence set pkey={pkey}"
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

def getSaltLetters(saSmileFragments):
    saSaltLetters = ''
    saRemainderSmailes = list()
    for smile in saSmileFragments:
        mol = Chem.MolFromSmiles(smile)
        # Remove all stereochemistry from the fragment
        Chem.rdmolops.RemoveStereochemistry(mol)
        flattenSmile = Chem.rdmolfiles.MolToSmiles(mol)
        canonSmile = Chem.CanonSmiles(flattenSmile)
        sSql = f"select suffix from chem_reg.salts where smiles='{flattenSmile}'"
        cur.execute(sSql)
        suffix = cur.fetchall()
        if len(suffix) == 0:
            saRemainderSmailes.append(smile)
        else:
            saSaltLetters += suffix[0][0]
    
    sRemainderSmile = '.'.join(saRemainderSmailes)
    if saSaltLetters == '' or len(saRemainderSmailes) > 1:
        # This is an error state, there are multiple fragments in
        # the molfile but there is no match against the salt database
        logger.error(f"Can't find salt in fragments {saSmileFragments}")
        return False, saRemainderSmailes

    return saSaltLetters, sRemainderSmile

def addStructure(database, molfile, newRegno, idColumnName):
    #####
    # Add the molecule to the structure tables
    sSql = f"""
    insert into {database} (mol, {idColumnName})
    value
    (mol2bin('{molfile}', 'mol'), '{newRegno}')
    """
    cur.execute(sSql)

    sSql = f"""
    insert into {database}_ukey select molid, uniquekey(mol) as molkey
    from {database} where {idColumnName} = '{newRegno}'
    """
    cur.execute(sSql)

    sSql = f"""
    insert into {database}_key select molid, fp(mol, 'sss') as molkey
    from {database} where {idColumnName} = '{newRegno}'
    """
    cur.execute(sSql)
    #
    #####

def registerNewCompound(compound_id_numeric,
                        molfile,
                        mf,
                        sep_mol_monoiso_mass,
                        ip_rights = '',
                        compound_name = ''):
    compound_id = f'CBK{compound_id_numeric}'
    sSql = f'''
    insert into bcpvs.compound (
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

def registerNewBatch(compound_id,
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
    sSql = f'''insert into bcpvs.batch (
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
    cur.execute(sSql)


def getMoleculeProperties(self, molfile):
    sio = sys.stderr = StringIO()
    Chem.WrapLogs()
    mol = Chem.MolFromMolBlock(molfile)
    try:
        C_MF = rdMolDescriptors.CalcMolFormula(mol)
        molmassFormula = Formula(C_MF.replace('-', ''))
        C_CHNS = getAtomicComposition(molmassFormula.composition())
    except Exception as e:
        return (False, False, False, False, False, '', f'{sio.getvalue()}')
    sSmiles = Chem.MolToSmiles(mol)
    if sSmiles == '':
        return (False, False, False, False, False, '', 'Empty molfile')
    saRemainderSmile = ''
    saSmileFragments = sSmiles.split('.')
    saSalts = ''
    if len(saSmileFragments) > 1:
        saSalts, saRemainderSmile = getSaltLetters(saSmileFragments)
        if saSalts == False:
            return (False,
                    False,
                    False,
                    False,
                    False,
                    '',
                    f'Unknown salt {saRemainderSmile}')
        #for smileFrag in saSmileFragments:
    if saRemainderSmile != '':
        resSmiles = saRemainderSmile
    else:
        resSmiles = sSmiles
    #######################
    ## Old code here
    C_MW = Descriptors.MolWt(mol)
    mono_iso_mol = Chem.MolFromSmiles(resSmiles)
    C_MONOISO = Descriptors.ExactMolWt(mono_iso_mol)
    return (C_MF, C_MW, C_MONOISO, C_CHNS, saSalts, resSmiles, '')

@jwtauth
class CreateSalt(tornado.web.RequestHandler):
    def put(self):
        sSmiles = self.get_argument("smiles")
        iSaltNumber = getNewSaltNumber()
        mol = Chem.MolFromSmiles(sSmiles)
        mw = Descriptors.MolWt(mol)
        mf = rdMolDescriptors.CalcMolFormula(mol)
        suffix = f'X{iSaltNumber}'
        sSql = f"""insert into chem_reg.salts (suffix, smiles, mw, mf)
        values ('{suffix}', '{sSmiles}', {mw}, '{mf}')
        """
        cur.execute(sSql)


@jwtauth
class GetCanonicSmiles(tornado.web.RequestHandler):
    def get(self):
        sSmiles = self.get_argument("smiles")
        sLetters, saRemainderSmile = getSaltLetters([sSmiles])
        mol = Chem.MolFromSmiles(sSmiles)
        flattenSmile = Chem.rdmolfiles.MolToSmiles(mol)
        canonSmile = Chem.CanonSmiles(flattenSmile)
        sRes = json.dumps([f'{sLetters}', canonSmile])
        self.finish(sRes)
        

@jwtauth
class GetLastBatchFromEln(tornado.web.RequestHandler):
    def get(self):
        sEln = self.get_argument("eln")
        eln_id = sEln[0:5]
        sSql = f'''select notebook_ref from bcpvs.batch
        where notebook_ref like '{eln_id}%\' order by notebook_ref desc'''
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
        sSql = "select pkey from chem_reg.sdfile_sequence"
        cur.execute(sSql)
        pkey = cur.fetchall()[0][0] +1
        sSql = f"update chem_reg.sdfile_sequence set pkey={pkey}"
        cur.execute(sSql)
        self.finish(f'{pkey}')


@jwtauth
class GetRegnosFromSdfSequence(tornado.web.RequestHandler):
    def get(self):
        sdfile_sequence = self.get_argument("sdfile_sequence")
        sSql = f'''select regno from chem_reg.chem_info
        where sdfile_sequence = {sdfile_sequence}
        '''
        res = sqlExec(sSql)
        self.finish(res)
        

@jwtauth
class BcpvsRegCompound(tornado.web.RequestHandler):
    def put(self):
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
        from chem_reg.chem_info
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
         sSmiles,
         errorMessage) = getMoleculeProperties(self, molfile)

        #mol = Chem.MolFromMolBlock(molfile)
        #sSmiles = Chem.MolToSmiles(mol)
        
        if compound_id in ('', None):
            mols = checkUniqueStructure(sSmiles)
            if len(mols) != 0:
                compound_id = mols[0][0]
                compound_id_numeric = compound_id[3:]
            else:
                # Register a new compound
                compound_id_numeric = getNewCompoundId()
                compound_id = registerNewCompound(compound_id_numeric,
                                                  molfile,
                                                  C_MF,
                                                  C_MONOISO,
                                                  ip_rights,
                                                  compound_name = '')
                addStructure("bcpvs.jcmol_moltable",
                             molfile,
                             compound_id,
                             'compound_id')
        sSql = f'''update chem_reg.chem_info
                   set compound_id='{compound_id}'
                   where regno = {regno}'''
        cur.execute(sSql)
        registerNewBatch(compound_id,
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
        molfile = self.get_body_argument('molfile')
        jpage = self.get_body_argument('jpage')
        chemist = self.get_body_argument('chemist')
        compound_type = self.get_body_argument('compound_type')
        project = self.get_body_argument('project')
        source = self.get_body_argument('source')
        solvent = self.get_body_argument('solvent')
        product = self.get_body_argument('product')
        library_id = self.get_body_argument('library_id')
        external_id = self.get_body_argument('external_id')
        supplier_batch = self.get_body_argument('supplier_batch')
        purity = self.get_body_argument('purity')
        ip_rights = self.get_body_argument('ip_rights')
        sdfile_sequence = self.get_body_argument('sdfile_sequence')
        
        (C_MF,
         C_MW,
         C_MONOISO,
         C_CHNS,
         saSalts,
         sSmiles,
         errorMessage) = getMoleculeProperties(self, molfile)

        if C_MF == False:
            self.set_status(500)
            self.finish(f'Molfile failed {external_id} {errorMessage} {sSmiles}')
            logger.error(f'Molfile failed {external_id} {errorMessage} {sSmiles}')
            return
        if saSalts == False:
            self.set_status(500)
            self.finish(f'Unknown salt in molfile for {external_id} {errorMessage} {sSmiles}')
            logger.error(f'Unknown salt in molfile for {external_id} {errorMessage} {sSmiles}')
            return

        ########################
        # Do exact match with molecule against the CBK database
        mols = checkUniqueStructure(sSmiles)
        sStatus = 'oldMolecule'
        if len(mols) == 0:
            sStatus = 'newMolecule'
            compound_id = ''
        else:
            compound_id = mols[0][0]

        ########################
        # Get new regno and add nmolecule to chem_info
        newRegno = getNewRegno()
        if purity == '':
            purity = -1
        sSql = f"""
        insert into chem_reg.chem_info (
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
        addStructure("chem_reg.chem_info_mol", molfile, newRegno, 'regno')
        self.finish(sStatus)
          

@jwtauth
class Search(tornado.web.RequestHandler):
    def get(self):
        column = self.get_argument("column")
        value = self.get_argument("value")
        values = (value, )

        sSql = "select regno from chem_reg.chem_info where " + column +" = %s"
        res = sqlExec(sSql, values)
        self.write(res)


@jwtauth
class GetRegnoData(tornado.web.RequestHandler):
    def get(self):
        column = self.get_argument("column")
        regno = self.get_argument("regno")
        values = (regno, )        
        sSql = "select " + column + " from chem_reg.chem_info where regno = %s"
        res = sqlExec(sSql, values)
        self.write(res)


@jwtauth
class CreateMolImage(tornado.web.RequestHandler):
    def get(self):
        regno = self.get_argument("regno")
        sSql = f"""select molfile from chem_reg.chem_info
                   where regno = '{regno}'"""
        cur.execute(sSql)
        molfile = cur.fetchall()
        if len(molfile) > 0 and molfile[0][0] != None:
            createPngFromMolfile(regno, molfile[0][0])
        self.finish()


@jwtauth
class LoadMolfile(tornado.web.RequestHandler):
    def post(self):
        fBody = self.request.files['file'][0]
        regnoBody = self.request.files['regno'][0]
        molfile = tornado.escape.xhtml_unescape(fBody.body)
        regno = tornado.escape.xhtml_unescape(regnoBody.body)
        (C_MF,
         C_MW,
         C_MONOISO,
         C_CHNS,
         saSalts,
         sSmiles,
         errorMessage) = getMoleculeProperties(self, molfile)
        if C_MF == False:
            self.set_status(500)
            self.finish(f'Molfile failed {regno} {errorMessage}')
            print(f'Molfile failed {regno}')
            return
        try:
            addStructure('chem_reg.chem_info_mol', molfile, regno, 'regno')
        except Exception as e:
            print(str(e))
            print(f'failed on regno {regno}')
            return
        sSql = f"""update chem_reg.chem_info set
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
        regno = self.get_argument("regno")
        batch = self.get_argument("batch")
        sSql = f"""select compound_id from bcpvs.batch
                   where notebook_ref = '{batch}'"""
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            self.set_status(500)
            self.finish()
        else:
            sSql = f"""select regno from chem_reg.chem_info
            where jpage = '{batch}' and regno != '{regno}'"""
            cur.execute(sSql)
            res = cur.fetchall()
            if len(res) > 0:
                self.set_status(500)
                self.finish()
            else:
                sSql = f"""update chem_reg.chem_info
                set jpage = '{batch}' where regno = '{regno}'"""
                cur.execute(sSql)
                self.finish()


@jwtauth
class CreateSupplier(tornado.web.RequestHandler):
    def put(self):
        supplierName = self.get_argument("supplier")
        sSql = f"""insert into bcpvs.compound_suppliers
        (name) values ('{supplierName}')"""
        cur.execute(sSql)


@jwtauth
class CreateRegno(tornado.web.RequestHandler):
    def put(self):
        regno = self.get_argument("regno")
        sSql = """insert into chem_reg.chem_info
        (regno, rdate) values (%s, now())"""
        val = (regno, )
        cur.execute(sSql, val)


@jwtauth
class DeleteRegno(tornado.web.RequestHandler):
    def put(self):
        regno = self.get_argument("regno")
        val = (regno, )
        sSql = """select regno, c_mf, chemist from chem_reg.chem_info
                  where c_mf is null and chemist is null and
                  regno=%s"""
        cur.execute(sSql, val)
        res = cur.fetchall()
        if len(res) > 0:
            logger.info('Deleting ' + str(res))
            sSql = """delete from chem_reg.chem_info
                      where regno = %s"""
            cur.execute(sSql, val)


@jwtauth
class UpdateColumn(tornado.web.RequestHandler):
    def put(self):
        column = self.get_argument("column")
        value = self.get_argument("value")
        regno = self.get_argument("regno")
        values = (value, regno, )
        sSql = "update chem_reg.chem_info set " + column
        sSql += """= %s where regno = %s"""
        cur.execute(sSql, values)


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
        regno = self.get_argument("regno")
        sSql = f"""select molfile from chem_reg.chem_info
                   where regno = '{regno}'"""
        cur.execute(sSql)
        res = cur.fetchall()
        if len(res) > 0:
            self.write(res[0][0])


@jwtauth
class GetTextColumn(tornado.web.RequestHandler):
    def get(self):
        column = self.get_argument("column")
        regno = self.get_argument("regno")
        sSql = "select cast(" + column + """ as char)
                from chem_reg.chem_info where regno = %s"""
        values = (regno, )
        res = sqlExec(sSql, values)
        self.write(res)


@jwtauth
class GetColComboData(tornado.web.RequestHandler):
    def get(self):
        column = self.get_argument("column")
        sSql = ''
        if column == 'project':
            sSql = """select project_name from hive.project_details
                      order by project_name"""
        elif column == 'chemist':
            sSql = """select fullname from hive.user_details
            where ORGANIZATION = 'chemistry'"""
        elif column == 'compound_type':
            sSql = "select type from bcpvs.compound_type order by type"
        elif column == 'product':
            sSql = "SELECT type FROM bcpvs.product_type order by type"
        elif column == 'supplier':
            sSql = "SELECT name FROM bcpvs.compound_suppliers order by name"
        elif column == 'solvent':
            sSql = "SELECT solvent FROM chemspec.solvent_tbl order by solvent"
        elif column == 'library_id':
            sSql = """SELECT library_name FROM bcpvs.compound_library
                      order by library_name"""
        elif column == 'library_description':
            sSql = """SELECT description FROM bcpvs.compound_library
                      order by description"""

        res = sqlExec(sSql)
        self.write(res)


@jwtauth
class GetLibraryName(tornado.web.RequestHandler):
    def get(self):
        library_id = self.get_argument("library_id")
        values = (library_id, )
        sSql = """SELECT description FROM bcpvs.compound_library
                  where library_name = %s"""
        res = sqlExec(sSql, values)
        self.write(res)


@jwtauth
class GetLibraries(tornado.web.RequestHandler):
    def get(self):
        sSql = "select fullname from hive.user_details where ORGANIZATION = 'chemistry'"
        res = sqlExec(sSql)
        self.write(res)


@jwtauth
class GetNextRegno(tornado.web.RequestHandler):
    def get(self):
        newRegno = getNewRegno()
        self.write(json.dumps(newRegno))


@jwtauth
class CreateLibrary(tornado.web.RequestHandler):
    def put(self):
        sLibraryDescription = self.get_argument("library_name")
        sSupplier = self.get_argument("supplier")

        sSql = "select pkey from bcpvs.library_id_sequence"
        cur.execute(sSql)
        pkey = cur.fetchall()[0][0] +1
        sSql = f"update bcpvs.library_id_sequence set pkey={pkey}"
        cur.execute(sSql)

        library_id = f'Lib-{pkey}'
        
        sSql = f"""insert into bcpvs.compound_library
        (library_name,
        supplier,
        description) values
        ('{library_id}',
        '{sSupplier}',
        '{sLibraryDescription}')
        """
        cur.execute(sSql)
