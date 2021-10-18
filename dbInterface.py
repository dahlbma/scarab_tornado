import tornado.web
import json
#import jwt
import MySQLdb
from auth import jwtauth
import config
import logging

logger = logging.getLogger(__name__)

db_connection = MySQLdb.connect(
    host=config.database['host'],
    user=config.database['user'],
    passwd=config.database['password']
)
db_connection.autocommit(True)
cur = db_connection.cursor()

def sqlExec(sSql):
    cur.execute(sSql)
    res = [list(i) for i in cur.fetchall()]
    return json.dumps(res)

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
        
@jwtauth
class CreateRegno(tornado.web.RequestHandler):
    def put(self):
        regno = self.get_argument("regno")
        sSql = """insert into chem_reg.chem_info (regno, rdate) values (%s, now())"""
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
