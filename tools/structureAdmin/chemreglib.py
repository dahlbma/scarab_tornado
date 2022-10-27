import sys, requests, json, dbInterface, os, subprocess, platform, shutil, datetime, traceback, logging
from PyQt5.QtGui import QImage, QPixmap
from PyQt5 import QtWidgets, QtCore
from PyQt5.QtWidgets import QApplication, QMessageBox

ip_rights_list = [' ', 'External rights', 'LCBKI', 'Commercial']


def send_msg(title, text, icon=QMessageBox.Information, e=None):
    msg = QMessageBox()
    msg.setWindowTitle(title)
    msg.setIcon(icon)
    msg.setText(text)
    clipboard = QApplication.clipboard()
    if e is not None:
        # add clipboard btn
        msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Save)
        buttonS = msg.button(QMessageBox.Save)
        buttonS.setText('Save to clipboard')
    msg.exec_()
    if e is not None:
        if msg.clickedButton() == buttonS:
            # copy to clipboard if clipboard button was clicked
            clipboard.setText(text)
            cb_msg = QMessageBox()
            cb_msg.setText(clipboard.text()+" \n\ncopied to clipboard!")
            cb_msg.exec_()

def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

def open_file(filename):
    # open file with default OS application
    # opens folders too
    if platform.system() == 'Windows':
        if filename == ".":
            proc = os.startfile(filename)
        else:
            proc = subprocess.Popen("start " + filename, shell=True)
    elif platform.system() == 'Linux':
        proc = subprocess.Popen("xdg-open " + filename, shell=True)
    elif platform.system() == 'Darwin':
        proc = subprocess.Popen("open " + filename, shell=True)
    return proc

def displayMolfile(self, sId):
    #dbInterface.createMolImage(self.token,
    #                           self.compound_id)
    image = QImage()
    self.new_structure_lab.setScaledContents(True)
    image.loadFromData(dbInterface.getMolImage(sId))
    self.new_structure_lab.setPixmap(QPixmap(image))

def postMolFile(self, fname, regno, logger):
    logger.info("posting file %s to server", fname)
    f = {'file': open(fname, 'rb'), 'regno': regno}
    r = dbInterface.postMolfile(self.token,
                                f)
    if r.status_code != 200:
        return False, r.content
    else:
        return True, r.content
    
