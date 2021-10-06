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

# Secret stuff in config file
import config

enable_pretty_logging()

db_connection = MySQLdb.connect(
    host=config.database['host'],
    user=config.database['user'],
    passwd=config.database['password']
)

cur = db_connection.cursor()
cur.execute( "SELECT * FROM chemspec.chem_info limit 1")
res = cur.fetchall()
print(res)

root = os.path.dirname(__file__)
settings = {
    "static_path": os.path.join(os.path.dirname(__file__), "angular-tour-of-heroes"),
    "cookie_secret": "__TODO:_GENERATE_YOUR_OWN_RANDOM_VALUE_HERE__",
}

JWT_SECRET = 'secret'
JWT_ALGORITHM = 'HS256'
JWT_EXP_DELTA_SECONDS = 20

data = [
    { "id":11, "name":"Dr Nice" },
    { "id":12, "name":"Narco" },
    { "id":13, "name":"Bombasto" },
    { "id":14, "name":"Celeritas" },
    { "id":15, "name":"Magneta" },
    { "id":16, "name":"RubberMan" },
    { "id":17, "name":"Dynama" },
    { "id":18, "name":"Dr IQ" },
    { "id":19, "name":"Magma" },
    { "id":20, "name":"Tornado" }
]

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
    '''
    def set_default_headers(self):
        self.set_header("Access-Control-Allow-Origin", "*")
        self.set_header("Access-Control-Allow-Headers", "x-requested-with")
        self.set_header('Access-Control-Allow-Methods', "GET,HEAD,OPTIONS,POST,PUT")
        self.set_header("Access-Control-Allow-Credentials", "true")
        self.set_header("Access-Control-Allow-Headers", "Access-Control-Allow-Headers,\
        Origin,Accept, X-Requested-With, Content-Type, Access-Control-Request-Method,\
        Access-Control-Request-Headers")
    '''
    def set_default_headers(self):
        self.set_header("Access-Control-Allow-Origin", "*")
        self.set_header("Access-Control-Allow-Headers", "x-requested-with")
        self.set_header('Access-Control-Allow-Methods', "GET,HEAD,OPTIONS,POST,PUT")
        self.set_header("Access-Control-Allow-Credentials", "true")
        self.set_header("Access-Control-Allow-Headers", "*")

    def post(self):
        self.write('some post')

    def get(self):
        self.write('some get')

    def options(self, *args):
        # no body
        # `*args` is for route with `path arguments` supports
        self.set_status(204)
        self.finish()

#class MainHandler(tornado.web.RequestHandler):
class MainHandler(BaseHandler):
    def get(self):
        #res = json.loads(data)
        self.write(json.dumps(data))

#class getHero(tornado.web.RequestHandler):
class getHero(BaseHandler):
    def get(self, hero_id):
        for i in data:
            if i['id'] == int(hero_id):
                print(json.dumps(i))
                self.write(json.dumps(i))

class getLogin(BaseHandler): 
    def post(self):
        error, username, password = getArgs(self.request.body)

        print(username)
        print(password)
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
        (r"/api/auth/signin", getLogin),
        (r"/api/heroes/([0-9]+)", getHero),
        (r"/heroes", MainHandler),
        (r"/", MainHandler),
        (r"/(.*)", tornado.web.StaticFileHandler, {"path": root,
                                                   "default_filename": "index.html"}),
    ], **settings)


if __name__ == "__main__":
    app = make_app()
    app.listen(8082)
    tornado.autoreload.start()
    
    for dir, _, files in os.walk('static'):
        [tornado.autoreload.watch(dir + '/' + f) for f in files if not f.startswith('.')]

    tornado.ioloop.IOLoop.current().start()
