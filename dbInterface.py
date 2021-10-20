import tornado.web
import json
import MySQLdb
from auth import jwtauth
import config
import logging
from rdkit import Chem
from rdkit.Chem import Draw

logger = logging.getLogger(__name__)

db_connection = MySQLdb.connect(
    host=config.database['host'],
    user=config.database['user'],
    passwd=config.database['password']
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
        fileName = "mols/" + regno + ".mol"
        fileHandle = open(fileName, "w")
        fileHandle.write(molfile)
        fileHandle.close()
        m = Chem.MolFromMolFile(fileName)
        Draw.MolToFile(m,'mols/' + regno + '.png')

        
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
        #print('Deleting')
        #print(res)
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
        sSql = "select " + column + " from chem_reg.chem_info where regno = %s"
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
        elif column == "":
            pass
            
        res = sqlExec(sSql)
        self.write(res)


@jwtauth
class GetChemists(tornado.web.RequestHandler):
    def get(self):
        sSql = """select fullname from hive.user_details
                  where ORGANIZATION = 'chemistry'"""
        res = sqlExec(sSql)
        self.write(res)

@jwtauth
class GetProjects(tornado.web.RequestHandler):
    def get(self):
        #sSql = """select project_name from hive.project_details
        #          order by created_date desc"""
        sSql = """select project_name from hive.project_details
                  order by project_name"""
        res = sqlExec(sSql)
        self.write(res)

@jwtauth
class GetCompoundTypes(tornado.web.RequestHandler):
    def get(self):
        sSql = "select type from bcpvs.compound_type order by type"
        res = sqlExec(sSql)
        self.write(res)

@jwtauth
class GetProductTypes(tornado.web.RequestHandler):
    def get(self):
        sSql = "SELECT type FROM bcpvs.product_type order by type"
        res = sqlExec(sSql)
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
        sSql = "select id from chemspec.regno_sequence"
        cur.execute(sSql)
        id = cur.fetchall()[0][0] +1
        sSql = "update chemspec.regno_sequence set id=" + str(id)
        cur.execute(sSql)
        self.write(json.dumps(id))
