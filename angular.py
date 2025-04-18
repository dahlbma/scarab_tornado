import json
import tornado.ioloop
import tornado.web
from tornado.escape import json_decode
import os
import re
import MySQLdb
from datetime import datetime, timedelta
import jwt
from tornado import web
import tornado.autoreload
from tornado.log import enable_pretty_logging
import chembl_export
import logging
import dbInterface
from auth import jwtauth

# Secret stuff in config file
import config

logger = logging.getLogger(__name__)
enable_pretty_logging()


root = os.path.dirname(__file__)
settings = {
    "cookie_secret": config.secret_key,
    "static_path": os.path.join(os.path.dirname(__file__), "dist"),
    "static_url_prefix": "/dist/",
}

JWT_SECRET = config.secret_key
JWT_ALGORITHM = 'HS256'
JWT_EXP_DELTA_SECONDS = 99999

class BaseHandler(tornado.web.RequestHandler):
    def set_default_headers(self):
        pass
        
    def post(self):
        self.write('some post')

    def get(self):
        self.write('some get')

    def options(self, *args):
        # no body
        # `*args` is for route with `path arguments` supports
        self.set_status(204)
        self.finish()


class login(tornado.web.RequestHandler):
    def post(self, *args):
        username = self.get_argument('username')
        password = self.get_argument('password')
        database = self.get_argument('database')
        try:
            db_connection2 = MySQLdb.connect(
                host="esox3",
                user=username,
                passwd=password
            )
            db_connection2.close()
        except Exception as ex:
            logger.error(str(ex))
            self.set_status(400)
            self.write({'message': 'Wrong username password'})
            self.finish()
            return
        payload = {
            'username': username,
            'exp': datetime.utcnow() + timedelta(seconds=JWT_EXP_DELTA_SECONDS)
        }
        jwt_token = jwt.encode(payload, JWT_SECRET, JWT_ALGORITHM)
        logging.info(f'Login: {username} Database; {database}')
        self.write({'token': jwt_token,
                    'database': database})

    def get(self):
        pass

class getVersionData(tornado.web.RequestHandler):
    def post(self):
        pass

    def get(self):
        # send verdat

        # tentative
        #os_name = self.get_argument('os_name')
        ## query db for version
        #sSql = f'''select version
        #           from `chem_reg`.`chemreg_dist`
        #           where os = {os_name};
        #        '''
        #cur.execute(sSql)
        #res = cur.fetchall()
        #self.write({"version": f'{res[0][0]}'})

        try:
            with open('./ver.dat', 'r') as f:
                self.write(json.load(f))
                return
        except Exception as e:
            logger.error(str(e))
            self.set_status(500)
            self.write({'message': 'ver.dat not available'})

class getChemRegBin(tornado.web.RequestHandler):
    def post(self):
        pass

    def get(self):
        # send chemReg
        os_name = self.get_argument('os_name')
        # tentative
        
        #
        #if not (os_name == 'Windows' or os_name == 'Linux' or os_name == 'Darwin'):
        #    # unsupported OS
        #    self.set_status(500)
        #    self.write({'message': 'OS not supported'})
        #    return
        #try:
        #    sSql = f'''select program
        #           from chem_reg.chemreg_dist
        #           where os = {os_name};
        #        '''
        #    cur.execute(sSql)
        #    res = cur.fetchall()
        #    logger.info("sending bin file")
        #    self.set_status(200)
        #    self.write(res[0][0])
        #except Exception as e:
        #    logger.error(f"Did not send bin file, error: {str(e)}")

        bin_file = ""
        if os_name == 'Windows':
            bin_file = os_name + "/chemreg.exe"
        elif os_name == 'Linux':
            bin_file = os_name + "/chemreg"
        elif os_name == 'Darwin':
            bin_file = os_name + "/chemreg"
        else:
            # unsupported OS
            self.set_status(500)
            self.write({'message': 'OS not supported'})
            return
        try:
            with open(bin_file, 'rb') as f:
                logger.info("sending bin file")
                self.set_status(200)
                self.write(f.read())
        except Exception as e:
            logger.error(f"Did not send bin file, error: {str(e)}")


def make_app():
    return tornado.web.Application([
        (r"/login", login),
        (r"/getDatabase", dbInterface.GetDatabase),
        (r"/pingDB", dbInterface.PingDB),
        (r"/api/bcpvsRegCompound", dbInterface.BcpvsRegCompound),
        (r"/api/chemRegAddMol", dbInterface.ChemRegAddMol),
        (r"/api/search", dbInterface.Search),
        (r"/api/loadMolfile", dbInterface.LoadMolfile),
        (r"/api/update", dbInterface.UpdateColumn),
        (r"/api/getRegnoData", dbInterface.GetRegnoData),
        (r"/api/createRegno", dbInterface.CreateRegno),
        (r"/api/updateRegnoBatch", dbInterface.UpdateRegnoBatch),
        (r"/api/createSalt", dbInterface.CreateSalt),
        (r"/api/addNostructMol", dbInterface.AddNostructMol),
        (r"/api/createSupplier", dbInterface.CreateSupplier),
        (r"/api/createLibrary", dbInterface.CreateLibrary),
        (r"/api/deleteRegno", dbInterface.DeleteRegno),
        (r"/api/getColComboData", dbInterface.GetColComboData),
        (r"/api/getTextColumn", dbInterface.GetTextColumn),
        (r"/api/getCanonicSmiles", dbInterface.GetCanonicSmiles),
        (r"/api/getLibraryName", dbInterface.GetLibraryName),
        (r"/api/getLibraries", dbInterface.GetLibraries),
        (r"/api/getLastBatchFromEln", dbInterface.GetLastBatchFromEln),
        (r"/api/getNextRegno", dbInterface.GetNextRegno),
        (r"/api/getBackwardCompound", dbInterface.GetBackwardCompound),
        (r"/api/getBackwardRegno", dbInterface.GetBackwardRegno),
        (r"/api/getMolkeyct", dbInterface.GetMolkeyct),
        (r"/api/getCompoundDuplicates", dbInterface.GetCompoundDuplicates),
        (r"/api/getForwardCompound", dbInterface.GetForwardCompound),
        (r"/api/getForwardRegno", dbInterface.GetForwardRegno),
        (r"/api/getNextSdfSequence", dbInterface.GetNextSdfSequence),
        (r"/api/getRegnosFromSequence", dbInterface.GetRegnosFromSdfSequence),
        (r"/api/getMolfile", dbInterface.GetMolfile),
        (r"/api/getMolfileBcpvs", dbInterface.GetMolfileBcpvs),
        (r"/api/getRegnoFromCompound", dbInterface.GetRegnoFromCompound),
        (r"/api/getCompoundFromRegno", dbInterface.GetCompoundFromRegno),
        (r"/api/createMolImage", dbInterface.CreateMolImage),
        (r"/api/createBcpvsMolImage", dbInterface.CreateBcpvsMolImage),
        (r"/api/updateStructureAdmin", dbInterface.UpdateStructureAdmin),
        (r"/api/createMolImageFromMolfile", dbInterface.CreateMolImageFromMolfile),
        (r"/getCompound", dbInterface.GetCompound),
        (r"/mols/(.*)", web.StaticFileHandler, {"path": "mols/"}),
        (r"/getVersionData", getVersionData), # dbInterface.
        (r"/getChemRegBin", getChemRegBin), # dbInterface.
        (r"/uploadBinary", dbInterface.UploadBinary), # upload
        (r"/getChemRegBinary/(?P<os_name>[^\/]+)", dbInterface.GetChemRegBinary), # upload
        (r"/uploadVersionNo", dbInterface.UploadVersionNo),
        (r"/uploadLauncher", dbInterface.UploadLauncher), # upload
        (r"/getMolsoft/(.*)", web.StaticFileHandler, {"path": "dist/molsoft.zip"}),
        (r"/getChemRegLauncher/Windows/(.*)", web.StaticFileHandler, {"path": "dist/launchers/Windows/"}),
        (r"/getChemRegLauncher/Linux/(.*)", web.StaticFileHandler, {"path": "dist/launchers/Linux/"}),
        (r"/getChemRegLauncher/Darwin/(.*)", web.StaticFileHandler, {"path": "dist/launchers/Darwin/"}),
        (r"/chemblExport/(?P<sRIDX>[^\/]+)/(?P<project>[^\/]+)/(?P<ELN>[^\/]+)/(?P<sBatches>[^\/]+)", dbInterface.ChemblExport),
        (r"/chemblExport", dbInterface.ChemblExport),
        (r"/(.*)", web.StaticFileHandler,  {"path": "dist", "default_filename": "index.html"}),
    ], **settings)

if __name__ == "__main__":
    app = make_app()
    app.listen(8081)
    tornado.autoreload.start()
    
    for dir, _, files in os.walk('static'):
        [tornado.autoreload.watch(dir + '/' + f) \
         for f in files if not f.startswith('.')]

    tornado.ioloop.IOLoop.current().start()
