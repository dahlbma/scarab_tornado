[1mdiff --git a/angular.py b/angular.py[m
[1mindex 7f34d00..a22b3fd 100644[m
[1m--- a/angular.py[m
[1m+++ b/angular.py[m
[36m@@ -46,6 +46,7 @@[m [mclass BaseHandler(tornado.web.RequestHandler):[m
         self.set_status(204)[m
         self.finish()[m
 [m
[32m+[m
 class login(tornado.web.RequestHandler):[m
     def post(self, *args):[m
         username = self.get_argument('username')[m
[1mdiff --git a/dbInterface.py b/dbInterface.py[m
[1mindex a94174d..6310d00 100644[m
[1m--- a/dbInterface.py[m
[1m+++ b/dbInterface.py[m
[36m@@ -1254,6 +1254,8 @@[m [mclass GetColComboData(tornado.web.RequestHandler):[m
         if column == 'project':[m
             sSql = """select project_name from hive.project_details[m
                       order by project_name"""[m
[32m+[m[32m            sSql = """select project_name from hive.project_details[m
[32m+[m[32m                      order by created_date desc"""[m
         elif column == 'chemist':[m
             sSql = """select fullname from hive.user_details[m
             where ORGANIZATION = 'chemistry'"""[m
