import sys
from PyQt5.uic import loadUi
from PyQt5.QtGui import QIntValidator
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QDialog, QApplication, QWidget, QMainWindow, QFileDialog
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
        self.onlyInt = QIntValidator()

        self.regno = dbInterface.getNextRegno(self.token)
        self.regno_lab.setText(self.regno)
        dbInterface.createNewRegno(self.regno, self.token)
        
        submitters = dbInterface.getSubmitters(self.token)
        self.submitter_cb.addItems(submitters)
        self.submitter_cb.currentTextChanged.connect(
            lambda x: self.changeEvent(x, 'CHEMIST'))
        
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
            lambda x: self.changeEvent(x, 'product'))

        self.externalid_eb.editingFinished.connect(
            lambda: self.changeEvent(self.externalid_eb.text(), 'EXTERNAL_ID'))

        self.externalbatch_eb.editingFinished.connect(
            lambda: self.changeEvent(self.externalbatch_eb.text(), 'SUPPLIER_BATCH'))

        self.batch_eb.editingFinished.connect(
            lambda: self.changeEvent(self.batch_eb.text(), 'JPAGE'))

        self.chrom_text.textChanged.connect(
            lambda: self.changeEvent(self.chrom_text.toPlainText(), 'CHROM_TEXT'))
        self.chrompurity_eb.textChanged.connect(
            lambda: self.changeEvent(self.chrompurity_eb.text(), 'CHROM_PURITY'))
        self.chrompurity_eb.setValidator(self.onlyInt)
        
        self.nmr_text.textChanged.connect(
            lambda: self.changeEvent(self.nmr_text.toPlainText(), 'NMR_TEXT'))
        self.nmrpurity_eb.textChanged.connect(
            lambda: self.changeEvent(self.nmrpurity_eb.text(), 'NMR_PURITY'))
        self.nmrpurity_eb.setValidator(self.onlyInt)
        
        self.ms_text.textChanged.connect(
            lambda: self.changeEvent(self.ms_text.toPlainText(), 'MS_TEXT'))
        self.mspurity_eb.textChanged.connect(
            lambda: self.changeEvent(self.mspurity_eb.text(), 'MS_PURITY'))
        self.mspurity_eb.setValidator(self.onlyInt)
        
        self.comments_text.textChanged.connect(
            lambda: self.changeEvent(self.comments_text.toPlainText(), 'COMMENTS'))

        self.loadmol_btn.clicked.connect(self.getMolFile)

    def getMolFile(self):
        fname = QFileDialog.getOpenFileName(self, 'Open file', 
                                                '.', "Molfiles (*.mol)")
        print(fname[0])
        files = {'regno': self.regno, 'file': open(fname[0], 'rb')}
        #r = requests.post('http://esox3.scilifelab.se:8082/api/loadMolfile',
        #                  params={'column': target, 'value': value, 'regno': regno},
        #                  headers={'token': token}, files=files)
        r = requests.post('http://esox3.scilifelab.se:8082/api/loadMolfile',
                          headers={'token': self.token}, files=files)


        
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
