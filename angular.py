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
        self.write({'token': jwt_token})

    def get(self):
        pass

class getVersionData(tornado.web.RequestHandler):
    def post(self):
        pass

    def get(self):
        # send verdat
        try:
            with open('./ver.dat', 'r') as f:
                self.write(json.loads(f.read()))
                return
        except Exception as e:
            print(str(e))
            self.set_status(500)
            self.write({'message': 'ver.dat not available'})

class getChemRegBin(tornado.web.RequestHandler):
    def post(self):
        pass

    def get(self):
        # send chemReg
        os_name = self.get_argument('os_name')
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
        with open(bin_file, 'rb') as f:
            self.set_status(200)
            self.write(f.read())


def make_app():
    return tornado.web.Application([
        (r"/login", login),
        (r"/api/bcpvsRegCompound", dbInterface.BcpvsRegCompound),
        (r"/api/chemRegAddMol", dbInterface.ChemRegAddMol),
        (r"/api/search", dbInterface.Search),
        (r"/api/loadMolfile", dbInterface.LoadMolfile),
        (r"/api/update", dbInterface.UpdateColumn),
        (r"/api/getRegnoData", dbInterface.GetRegnoData),
        (r"/api/createRegno", dbInterface.CreateRegno),
        (r"/api/updateRegnoBatch", dbInterface.UpdateRegnoBatch),
        (r"/api/createSalt", dbInterface.CreateSalt),
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
        (r"/api/getNextSdfSequence", dbInterface.GetNextSdfSequence),
        (r"/api/getRegnosFromSequence", dbInterface.GetRegnosFromSdfSequence),
        (r"/api/getMolfile", dbInterface.GetMolfile),
        (r"/api/createMolImage", dbInterface.CreateMolImage),
        (r"/getCompound", dbInterface.GetCompound),
        (r"/mols/(.*)", web.StaticFileHandler, {"path": "mols/"}),
        (r"/getVersionData", getVersionData),
        (r"/getChemRegBin", getChemRegBin)
    ], **settings)

if __name__ == "__main__":
    app = make_app()
    app.listen(8082)
    tornado.autoreload.start()
    
    for dir, _, files in os.walk('static'):
        [tornado.autoreload.watch(dir + '/' + f) \
         for f in files if not f.startswith('.')]

    tornado.ioloop.IOLoop.current().start()
