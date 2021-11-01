import tornado.web
import json
import MySQLdb
from auth import jwtauth
import config
import logging
from rdkit import Chem
from rdkit.Chem import Draw
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

def createPngFromMolfile(regno, molfile):
    fileName = "mols/" + regno + ".mol"
    fileHandle = open(fileName, "w")
    fileHandle.write(molfile)
    fileHandle.close()
    m = Chem.MolFromMolFile(fileName)
    try:
        Draw.MolToFile(m,'mols/' + regno + '.png')
    except:
        '''
        # Check:
        https://sourceforge.net/p/rdkit/mailman/rdkit-discuss/thread/48DA553F-AA94-455E-98D7-71287037FE6F%40gmail.com/
        '''
        logger.info(f"regno {regno} is nostruct")


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

        ####
        # Get pkey for tmp_mol table
        sSql = "select pkey from chem_reg.tmp_mol_sequence"
        cur.execute(sSql)
        pkey = cur.fetchall()[0][0] +1
        sSql = f"update chem_reg.tmp_mol_sequence set pkey={pkey}"
        cur.execute(sSql)

        def to_bytes(s):
            if type(s) is bytes:
                return s
            elif type(s) is str or type(s) is unicode:
                return codecs.encode(s, 'utf-8')
            else:
                raise TypeError("Expected bytes or string, but got %s." % type(s))

        ####
        # Insert molfile in tmp_mol table
        sSql = f"""
        insert into chem_reg.tmp_mol (pkey, molfile) values
        ({pkey}, '{molfile}')
        """
        cur.execute(sSql)
        
        ####
        # Do exact match with molecule against present molucules
        sSql = f"""
        select bin2smiles(chem_reg.mol.mol) from
          chem_reg.mol_ukey join mol on (chem_reg.mol.molid=chem_reg.mol_ukey.molid)
        where uniquekey(mol2bin(
            'select molfile from chem_reg.tmp_mol where pkey={pkey}', 'mol'))=molkey
        """
        cur.execute(sSql)
        mols = cur.fetchall()
        
        ####
        # Get new regno
        newRegno = getNewRegno()
        if purity == '':
            purity = -1
        sSql = f"""
        insert into chem_reg.chem_info (
        regno,
        jpage,
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
        molfile)
        values (
        '{newRegno}',
        '{jpage}',
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
        '{molfile}'
        )
        """
        cur.execute(sSql)
        createPngFromMolfile(str(newRegno), molfile)
        
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
            
        ####
        # Cleanup tmp_mol table, delete the temporary molfile
        sSql = f"""delete from chem_reg.tmp_mol where pkey={pkey}"""
        cur.execute(sSql)


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
class LoadMolfile(tornado.web.RequestHandler):
    def post(self):
        fBody = self.request.files['file'][0]
        regnoBody = self.request.files['regno'][0]
        molfile = tornado.escape.xhtml_unescape(fBody.body)
        regno = tornado.escape.xhtml_unescape(regnoBody.body)
        sSql = """update chem_reg.chem_info set `molfile` = %s
                  where regno = %s"""
        values = (molfile, regno, )
        cur.execute(sSql, values)
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
        logger.info('Deleting ' + str(res))
        if len(res) > 0:
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
