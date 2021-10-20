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

def getArgs(sBody):
    error = False
    username = ''
    password = ''
    
    data = tornado.escape.xhtml_unescape(sBody)
    x = re.match("^{username=(.+)&password=(.+)}$", data)

    if len(x.groups()) != 2:
        error = True
        
    if error != True:
        username = x.groups()[0]
        password = x.groups()[1]
    return (error, username, password)

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
        
class getLogin(BaseHandler):
    def post(self):
        error, username, password = getArgs(self.request.body)
        try:
            db_connection2 = MySQLdb.connect(
                host="esox3",
                user=username,
                passwd=password
            )
        except:
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
        
def make_app():
    return tornado.web.Application([
        (r"/login", login),
        (r"/api/search", dbInterface.Search),
        (r"/api/loadMolfile", dbInterface.LoadMolfile),
        (r"/api/update", dbInterface.UpdateColumn),
        (r"/api/getRegnoData", dbInterface.GetRegnoData),
        (r"/api/createRegno", dbInterface.CreateRegno),
        (r"/api/deleteRegno", dbInterface.DeleteRegno),
        (r"/api/getColComboData", dbInterface.GetColComboData),
        (r"/api/getTextColumn", dbInterface.GetTextColumn),
        (r"/api/getChemists", dbInterface.GetChemists),
        (r"/api/getProjects", dbInterface.GetProjects),
        (r"/api/getCompoundTypes", dbInterface.GetCompoundTypes),
        (r"/api/getProductTypes", dbInterface.GetProductTypes),
        (r"/api/getLibraries", dbInterface.GetLibraries),
        (r"/api/getNextRegno", dbInterface.GetNextRegno),
        (r"/getCompound", dbInterface.GetCompound),
        (r"/api/auth/signin", getLogin),
        (r"/mols/(.*)", web.StaticFileHandler, {"path": "mols/"}),
    ], **settings)

if __name__ == "__main__":
    app = make_app()
    app.listen(8082)
    tornado.autoreload.start()
    
    for dir, _, files in os.walk('static'):
        [tornado.autoreload.watch(dir + '/' + f) \
         for f in files if not f.startswith('.')]

    tornado.ioloop.IOLoop.current().start()
