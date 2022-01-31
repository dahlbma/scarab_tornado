import sys, requests, logging
from PyQt5.uic import loadUi
from PyQt5 import QtWidgets, QtCore
from PyQt5.QtWidgets import QDialog
from PyQt5.QtWidgets import QMessageBox

from chemreglib import *
from regscreen import RegScreen
from searchscreen import SearchScreen


class LoginScreen(QDialog):
    def __init__(self):
        super(LoginScreen, self).__init__()
        self.mod_name = "login"
        logger = logging.getLogger(self.mod_name)
        loadUi(resource_path("assets/welcomescreen.ui"), self)
        self.passwordfield.setEchoMode(QtWidgets.QLineEdit.Password)
        self.login.clicked.connect(self.loginfunction)
        saDatabases = dbInterface.getDatabase()
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
        self.gotoSearch(self.jwt_token)

    def gotoReg(self, token):
        reg = RegScreen(token)
        self.window().addWidget(reg)
        self.window().setCurrentIndex(self.window().currentIndex() + 1)

    def gotoSearch(self, token):
        search = SearchScreen(token)
        self.window().addWidget(search)
        self.window().setCurrentIndex(self.window().currentIndex() + 1)

