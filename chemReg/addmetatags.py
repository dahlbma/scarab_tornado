import dbInterface, os, datetime, re, codecs, logging
from PyQt5.uic import loadUi
from PyQt5 import QtCore
from PyQt5.QtWidgets import QDialog, QApplication
from PyQt5.QtWidgets import QFileDialog, QProgressBar, QMessageBox

from chemreglib import *

class AddMetaTags(QDialog):
    def __init__(self, token):
        super(AddMetaTags, self).__init__()
        self.token = token
        self.mod_name = "addmetatags"
        logger = logging.getLogger(self.mod_name)
        loadUi(resource_path("assets/addmetatags.ui"), self)
        suppliers = dbInterface.getColComboData(self.token, 'supplier')
        self.supplier_cb.addItems(suppliers)
        self.supplier_cb.setCurrentText(' ')
        
        self.wipe = True
        self.tabWidget.setCurrentIndex(0)
        self.tabWidget.currentChanged.connect(self.tab_changed)
        
        #connect buttons
        self.library_save_btn.clicked.connect(self.library_save)
        self.supplier_save_btn.clicked.connect(self.supplier_save)
        self.salt_save_btn.clicked.connect(self.salt_save)
        
        self.exec_()
        self.show()

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Return or event.key() == QtCore.Qt.Key_Enter:
            i = self.tabWidget.currentIndex()
            if i == 0:
                self.library_save()
            elif i == 1:
                self.supplier_save()
            elif i == 2:
                self.salt_save()

    def library_save(self):
        lib_name = self.library_name_eb.text().strip()
        sup = self.supplier_cb.currentText()
        if lib_name == "":
            send_msg("Bad input", "Please enter a library name.")
            return
        # check if unique
        lib_low = lib_name.casefold()
        ll_names = dbInterface.getColComboData(self.token, 'library_description')
        for i in range(len(ll_names)):
            if ll_names[i] is None:
                continue
            ll_names[i] = ll_names[i].casefold()
        if lib_low in ll_names:
            send_msg("Library already exists", f"A library with the same name: \'{lib_name}\' already exists (case insensitive).")
            return
        
        ok = self.confirm_dialog("Do you want to add this library?" + f"\n Library name: {lib_name}" + f"\n Supplier: {sup}")
        if ok:
            #send
            check, r = dbInterface.createLibrary(self.token, lib_name, sup)
            if check:
                send_msg("Library added", f"Library: \'{lib_name}\', Supplier: \'{sup}\' added.")
                self.library_name_eb.setText(None)
                self.supplier_cb.setCurrentText(' ')
            else:
                logging.getLogger(self.mod_name).error(r.content)
                send_msg("Library NOT added", f"error: {r.content}")
    
    def supplier_save(self):
        sup_name = self.supplier_eb.text().strip()
        if sup_name == "":
            send_msg("Bad input", "Please enter a supplier name.")
            return
        
        # check if unique
        sup_low = sup_name.casefold()
        supps = dbInterface.getColComboData(self.token, 'supplier')
        for i in range(len(supps)):
            if supps[i] is None:
                continue
            supps[i] = supps[i].casefold()
        if sup_low in supps:
            send_msg("Supplier already exists", f"A supplier with the same name: \'{sup_name}\' already exists (case insensitive).")
            return
        
        ok = self.confirm_dialog("Do you want to add this supplier?" + f"\n Supplier: {sup_name}")
        if ok:
            #send
            check, r = dbInterface.createSupplier(self.token, sup_name)
            if check:
                send_msg("Supplier added", f"Supplier: \'{sup_name}\' added.")
                self.supplier_cb.addItem(sup_name)
                self.supplier_eb.setText(None)
            else:
                logging.getLogger(self.mod_name).error(r.content)
                send_msg("Supplier NOT added", f"error: {r.content}")
                
    
    def salt_save(self):
        salt = self.salt_eb.text().strip()
        if salt == "":
            send_msg("Bad input", "Please enter a salt.")
            return
        
        # check if unique
        canon_list = json.loads(dbInterface.getCanonicSmiles(self.token, salt))
        if canon_list[0] != 'False':
            send_msg("Salt already exists", f"A salt with the same canonic smiles: \'{canon_list[1]}\' already exists.")
        else:
            ok = self.confirm_dialog("Found a salt corresponding to input smiles-string. Do you want to add this salt?" + f"\n Salt: {canon_list[1]}")
            if ok:
                #send canonic smiles
                dbInterface.createSalt(self.token, canon_list[1])
                send_msg("Salt added", f"Salt: \'{canon_list[1]}\' added.")
                self.salt_eb.setText(None)
    
    def confirm_dialog(self, text):
        # confirm dialogue
        msg = QMessageBox()
        msg.setWindowTitle("Confirm Data")
        msg.setIcon(QMessageBox.Question)
        msg.setText(text)
        msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        msg.setDefaultButton(QMessageBox.Ok)
        btnS = msg.button(QMessageBox.Ok)
        msg.exec_()
        if msg.clickedButton() == btnS:
            return True
        else:
            return False
        
    def tab_changed(self):
        page_index = self.tabWidget.currentIndex()
        if self.wipe:
            if page_index == 0:
                self.supwipe()
                self.salwipe()
            elif page_index == 1:
                self.libwipe()
                self.salwipe()
            elif page_index == 2:
                self.libwipe()
                self.supwipe()

    def libwipe(self):
        self.library_name_eb.setText(None)
        self.supplier_cb.setCurrentText(' ')

    def supwipe(self):
        self.supplier_eb.setText(None)

    def salwipe(self):
        self.salt_eb.setText(None)
