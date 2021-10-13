import sys
from PyQt5.uic import loadUi
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QDialog, QApplication, QWidget, QMainWindow
import mysql
from mysql.connector import connect, Error
import requests

jwt_token = ''

class LoginScreen(QDialog):
    def __init__(self):
        super(LoginScreen, self).__init__()
        loadUi("welcomescreen.ui", self)
        self.passwordfield.setEchoMode(QtWidgets.QLineEdit.Password)
        self.login.clicked.connect(self.loginfunction)
        
    def loginfunction(self):
        user = self.usernamefield.text()
        password = self.passwordfield.text()

        if len(user) == 0 or len(password) == 0:
            self.errorlabel.setText("Please input all fields")
        else:
            self.errorlabel.setText("")

        r = requests.post('http://esox3.scilifelab.se:8082/login',
                          data = {'username':user, 'password':password})
        if r.status_code != 200:
            self.errorlabel.setText("Wrong username/password")
            return
        self.jwt_token = r.content
        self.gotoReg(self.jwt_token)

    def gotoReg(self, token):
        reg = RegScreen(token)
        widget.addWidget(reg)
        widget.setCurrentIndex(widget.currentIndex() + 1)

    def gotoSearch(self, token):
        search = SearchScreen()
        widget.addWidget(search)
        widget.setCurrentIndex(widget.currentIndex() + 1)


class RegScreen(QMainWindow):
    def __init__(self, token):
        super(RegScreen, self).__init__()
        loadUi("regchem.ui", self)
        self.token = token
        self.gotosearch_btn.clicked.connect(self.gotoSearch)
        self.setWindowTitle("Register new compound")
        print(self.token)
        r = requests.get('http://esox3.scilifelab.se:8082/api/getChemists',
                         headers={'token': self.token})

    def gotoSearch(self):
        search = SearchScreen()
        widget.addWidget(search)
        widget.setCurrentIndex(widget.currentIndex() + 1)
    

class SearchScreen(QMainWindow):
    def __init__(self):
        super(SearchScreen, self).__init__()
        loadUi("searchchem.ui", self)
        self.gotoreg_btn.clicked.connect(self.gotoReg)
        self.setWindowTitle("Search")
        self.regno_eb.setText('regno')

    def gotoReg(self):
        reg = RegScreen()
        widget.addWidget(reg)
        widget.setCurrentIndex(widget.currentIndex() + 1)

        
app = QApplication(sys.argv)
welcome = LoginScreen()
widget = QtWidgets.QStackedWidget()
widget.addWidget(welcome)
widget.setFixedHeight(800)
widget.setFixedWidth(1200)

widget.show()
try:
    sys.exit(app.exec_())
except:
    print("Exiting")
