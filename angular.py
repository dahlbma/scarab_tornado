import json
import tornado.ioloop
import tornado.web
from tornado.escape import json_decode
import os
import re
from urllib.parse import urlparse, parse_qs
from datetime import datetime, timedelta
import MySQLdb
import jwt
import tornado.autoreload
from tornado.log import enable_pretty_logging
from auth import jwtauth

# Secret stuff in config file
import config

enable_pretty_logging()

db_connection = MySQLdb.connect(
    host=config.database['host'],
    user=config.database['user'],
    passwd=config.database['password']
)

cur = db_connection.cursor()
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
    
        # self.set_header('Access-Control-Allow-Origin', '*')
        # self.set_header('Access-Control-Allow-Headers', '*')
        # self.set_header('Access-Control-Max-Age', 1000)
        # self.set_header('Content-type', 'application/json')
        # self.set_header('Access-Control-Allow-Methods', 'POST, GET, OPTIONS')
        # self.set_header('Access-Control-Allow-Headers', \
        #                 'Content-Type, Access-Control-Allow-Origin, \
        #                 Access-Control-Allow-Headers, X-Requested-By, \
        #                 Access-Control-Allow-Methods')
        
    def post(self):
        self.write('some post')

    def get(self):
        self.write('some get')

    def options(self, *args):
        # no body
        # `*args` is for route with `path arguments` supports
        self.set_status(204)
        self.finish()

@jwtauth
class RegCompound(tornado.web.RequestHandler):
    def get(self):
        # Contains user found in previous auth
        if self.request.headers.get('auth'):
            self.write('ok')

@jwtauth
class GetCompound(BaseHandler):
    def get(self):
        # Contains user found in previous auth
        if self.request.headers.get('auth'):
            self.write('ok')

    #get = post

#class MainHandler(tornado.web.RequestHandler):
class MainHandler(BaseHandler):
    def get(self):
        self.write(json.dumps(data))

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

def make_app():
    return tornado.web.Application([
        (r"/api/getCompound", GetCompound),
        (r"/api/auth/signin", getLogin),
        (r"/getCompound", GetCompound),
        (r"/(.*)", tornado.web.StaticFileHandler, {"path": "cbcs-compounds",
                                                   "default_filename": "index.html"}),
    ], **settings)

if __name__ == "__main__":
    app = make_app()
    app.listen(8082)
    tornado.autoreload.start()
    
    for dir, _, files in os.walk('static'):
        [tornado.autoreload.watch(dir + '/' + f) for f in files if not f.startswith('.')]

    tornado.ioloop.IOLoop.current().start()
