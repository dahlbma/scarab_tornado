import sys, requests, logging
from PyQt5.uic import loadUi
from PyQt5 import QtWidgets, QtCore
from PyQt5.QtWidgets import QDialog
from PyQt5.QtWidgets import QMessageBox

from chemreglib import *
from regscreen import RegScreen
from searchscreen import SearchScreen


class LoginScreen(QDialog):
    def __init__(self, appName):
        super(LoginScreen, self).__init__()
        self.mod_name = "login"
        self.appName = appName

        logger = logging.getLogger(self.mod_name)
        loadUi(resource_path("assets/welcomescreen.ui"), self)
        self.passwordfield.setEchoMode(QtWidgets.QLineEdit.Password)
        self.login.clicked.connect(self.loginfunction)
        saDatabases = None
        try:
            saDatabases = dbInterface.getDatabase()
        except Exception as e:
            send_msg("Connection Error", f"ChemReg has encountered a fatal error:\n\n{str(e)}\n\nPlease restart ChemReg.", icon=QMessageBox.Critical, e=e)
            sys.exit()
        self.server_cb.addItems(saDatabases)
        
    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Return \
           or event.key() == QtCore.Qt.Key_Enter:
            self.loginfunction()
    
    def loginfunction(self):
        user = self.usernamefield.text()
        password = self.passwordfield.text()
        database = self.server_cb.currentText()
        if len(user) == 0 or len(password) == 0:
            self.errorlabel.setText("Please input all fields")
        else:
            self.errorlabel.setText("")
        
        try:
            r = dbInterface.login(user, password, database)
        except Exception as e:
            self.errorlabel.setText("Bad Connection")
            send_msg("Error Message", str(e), QMessageBox.Warning, e)
            logging.getLogger(self.mod_name).error(str(e))
            return
        if r.status_code != 200:
            self.errorlabel.setText("Wrong username/password")
            return
        self.jwt_token = r.content
        app = QtCore.QCoreApplication.instance()
        self.window().setWindowTitle(f"{app.applicationName()} {database}")
        self.gotoSearch(self.jwt_token)

    def gotoReg(self, token):
        reg = RegScreen(token)
        self.window().addWidget(reg)
        self.window().setCurrentIndex(self.window().currentIndex() + 1)

    def gotoSearch(self, token):
        search = SearchScreen(token)
        self.window().setWindowTitle(f"{self.appName}")
        self.window().addWidget(search)
        self.window().setCurrentIndex(self.window().currentIndex() + 1)

