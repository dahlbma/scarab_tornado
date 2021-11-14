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
from rdkit.Chem.SaltRemover import SaltRemover
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


# Read the salt file
#with open('salts.json') as json_file:
#    salts = json.load(json_file)
sSql = f'select pkey, suffix, smiles, mf, mw from salts order by pkey'
cur.execute(sSql)
salts = cur.fetchall()

def sqlExec(sSql, values=None):
    if values == None:
        cur.execute(sSql)
    else:
        cur.execute(sSql, values)
    result = [list(i) for i in cur.fetchall()]
    return json.dumps(result)

def getNewRegno():
    sSql = "select id from chemspec.regno_sequence"
    cur.execute(sSql)
    id = cur.fetchall()[0][0] +1
    sSql = "update chemspec.regno_sequence set id=" + str(id)
    cur.execute(sSql)
    return id

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
        '''
        # Check:
        https://sourceforge.net/p/rdkit/mailman/rdkit-discuss/thread/48DA553F-AA94-455E-98D7-71287037FE6F%40gmail.com/
        '''
        logger.info(f"regno {regno} is nostruct")

def getSaltLetters(saSmileFragments):
    saSaltLetters = ''
    for smile in saSmileFragments:
        mol = Chem.MolFromSmiles(smile)
        # Remove all stereochemistry from the fragment
        Chem.rdmolops.RemoveStereochemistry(mol)
        flattenSmile = Chem.rdmolfiles.MolToSmiles(mol)
        canonSmile = Chem.CanonSmiles(flattenSmile)
        for salt in salts:
            saltMol = Chem.MolFromSmiles(salt[2])
            Chem.rdmolops.RemoveStereochemistry(saltMol)
            saltSmile = salt[2]
            if canonSmile == saltSmile:
                saSaltLetters += salt[1]
                break

    if saSaltLetters == '':
        # This is an error state, there are multiple fragments in
        # the molfile but there is no match against the salt database
        logger.info(f"Can't find salt in fragments {saSmileFragments}")
        return False
    return saSaltLetters

    
def getMoleculeProperties(self, molfile):
    sio = sys.stderr = StringIO()
    Chem.WrapLogs()
    mol = Chem.MolFromMolBlock(molfile)
    try:
        C_MF = rdMolDescriptors.CalcMolFormula(mol)
        molmassFormula = Formula(C_MF.replace('-', ''))
        C_CHNS = getAtomicComposition(molmassFormula.composition())
    except Exception as e:
        print(str(e))
        return (False, False, False, False, False, '', f'{sio.getvalue()}')
    sSmiles = Chem.MolToSmiles(mol)
    if sSmiles == '':
        return (False, False, False, False, False, '', 'Empty molfile')
    C_MW = Descriptors.MolWt(mol)
    C_MONOISO = Descriptors.ExactMolWt(mol)

    remover = SaltRemover()
    res = remover.StripMol(mol)
    numAtoms = mol.GetNumAtoms()
    salt = []
    saSalts = ''
    if res.GetNumAtoms() != numAtoms:
        res_mw = Descriptors.MolWt(res)
            
        saltMass = C_MW - res_mw
        sSmiles = Chem.MolToSmiles(mol)
        saSmileFragments = sSmiles.split('.')
        if len(saSmileFragments) > 1:
            saSalts = getSaltLetters(saSmileFragments)
    return (C_MF, C_MW, C_MONOISO, C_CHNS, saSalts, sSmiles, '')

        
@jwtauth
class chemRegAddMol(tornado.web.RequestHandler):
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
            logger.info(f'Molfile failed {external_id} {errorMessage} {sSmiles}')
            return
        if saSalts == False:
            self.set_status(500)
            self.finish(f'Unknown salt in molfile for {external_id} {errorMessage} {sSmiles}')
            logger.info(f'Unknown salt in molfile for {external_id} {errorMessage} {sSmiles}')
            return

        def to_bytes(s):
            if type(s) is bytes:
                return s
            elif type(s) is str or type(s) is unicode:
                return codecs.encode(s, 'utf-8')
            else:
                raise TypeError("Expected bytes or string, but got %s." % type(s))

        ####
        # Get new regno
        newRegno = getNewRegno()
        if purity == '':
            purity = -1
        sSql = f"""
        insert into chem_reg.chem_info (
        regno,
        jpage,
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
        C_CHNS,
        C_MF,
        C_MW,
        C_MONOISO,
        SUFFIX,
        molfile)
        values (
        '{newRegno}',
        '{jpage}',
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
        '{C_CHNS}',
        '{C_MF}',
        {C_MW},
        {C_MONOISO},
        '{saSalts}',
        '{molfile}'
        )
        """
        cur.execute(sSql)
        
        ####
        # Do exact match with molecule against present molucules
        sSql = f"""
        select bin2smiles(bcpvs.jcsepmol_moltable.mol) from
          bcpvs.jcsepmol_moltable_ukey join bcpvs.jcsepmol_moltable on
                 (bcpvs.jcsepmol_moltable.molid=bcpvs.jcsepmol_moltable_ukey.molid)
        where uniquekey(mol2bin('{molfile}', 'mol'))=molkey
        """
        cur.execute(sSql)
        mols = cur.fetchall()
        ####
        # Reg the molfile in chem_info if the molfile is unique        
        if len(mols) == 0:
            sSql = f"""
            insert into chem_reg.mol (mol, regno)
            value
            (mol2bin('{molfile}', 'mol'), {newRegno})
            """
            cur.execute(sSql)

            sSql = f"""
            insert into chem_reg.mol_ukey select molid, uniquekey(mol) as molkey
            from chem_reg.mol where regno = '{newRegno}'
            """
            cur.execute(sSql)

            sSql = f"""
            insert into chem_reg.mol_key select molid, fp(mol, 'sss') as molkey
            from chem_reg.mol where regno = '{newRegno}'
            """
            cur.execute(sSql)
            self.finish('newMolecule')
        else:
            self.finish('oldMolecule')
            

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
        sSql = f"""update chem_reg.chem_info set
                   `molfile` = '{molfile}',
                    C_MF = '{C_MF}',
                    C_MW = {C_MW},
                    C_MONOISO = {C_MONOISO},
                    C_CHNS = '{C_CHNS}'
                  where regno = {regno}"""
        cur.execute(sSql)
        createPngFromMolfile(regno, molfile)
        

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
