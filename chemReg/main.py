import sys
from PyQt5.uic import loadUi
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QDialog, QApplication, QWidget, QMainWindow
import mysql
from mysql.connector import connect, Error
import requests
import json
import dbInterface

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
        search = SearchScreen(token)
        widget.addWidget(search)
        widget.setCurrentIndex(widget.currentIndex() + 1)


class RegScreen(QMainWindow):
    def __init__(self, token):
        super(RegScreen, self).__init__()
        loadUi("regchem.ui", self)
        self.token = token
        self.dirty = False
        self.gotosearch_btn.clicked.connect(self.gotoSearch)
        self.setWindowTitle("Register new compound")

        self.regno = dbInterface.getNextRegno(self.token)
        self.regno_lab.setText(self.regno)
        dbInterface.createNewRegno(self.regno, self.token)
        
        submitters = dbInterface.getSubmitters(self.token)
        self.submitter_cb.addItems(submitters)
        self.submitter_cb.currentTextChanged.connect(
            lambda x: self.changeEvent(x, 'submitter'))
        
        projects = dbInterface.getProjects(self.token)
        self.project_cb.addItems(projects)
        self.project_cb.currentTextChanged.connect(
            lambda x: self.changeEvent(x, 'project'))
        
        compoundTypes = dbInterface.getCompoundTypes(self.token)
        self.compoundtype_cb.addItems(compoundTypes)
        self.compoundtype_cb.currentTextChanged.connect(
            lambda x: self.changeEvent(x, 'compound_type'))
        
        productTypes = dbInterface.getProductTypes(self.token)
        self.product_cb.addItems(productTypes)
        self.product_cb.currentTextChanged.connect(
            lambda x: self.changeEvent(x, 'produc'))
        
        #projects = dbInterface.getProjects(self.token)
        #self.project_cb.addItems(projects)
        #self.project_cb.currentTextChanged.connect(
        #    lambda x: self.changeEvent(x, 'updateProject'))

    def changeEvent(self, value='', action='doNothing'):
        if action == 'doNothing':
            return
        self.dirty = True
        dbInterface.updateValue(action, value, self.token, self.regno)

    def gotoSearch(self):
        if self.dirty == False:
            dbInterface.deleteRegno(self.regno, self.token)
        search = SearchScreen(self.token)
        widget.addWidget(search)
        widget.setCurrentIndex(widget.currentIndex() + 1)
    

class SearchScreen(QMainWindow):
    def __init__(self, token):
        super(SearchScreen, self).__init__()
        self.token = token
        loadUi("searchchem.ui", self)
        self.gotoreg_btn.clicked.connect(self.gotoReg)
        self.setWindowTitle("Search")
        self.regno_eb.setText('regno')

    def gotoReg(self):
        reg = RegScreen(self.token)
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
