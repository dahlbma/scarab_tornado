import sys, os
from PyQt5.uic import loadUi
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow

from chemreglib import *
from loginscreen import LoginScreen

def error_handler(etype, value, tb):
    err_msg = "".join(traceback.format_exception(etype, value, tb))
    #print("error caught")
    logging.getLogger().exception(err_msg + " from error handler\n")

logger = setLogger("mainLogger")
sys.excepthook = error_handler

try:
    os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "2"
    app = QApplication(['Chem Reg'])
    clipboard = app.clipboard()

    welcome = LoginScreen()
    widget = QtWidgets.QStackedWidget()
    widget.addWidget(welcome)

    desktop = QApplication.desktop()
    windowHeight = 800
    windowWidth = 1200

    windowHeight = int(round(0.9 * desktop.screenGeometry().height(), -1))
    if windowHeight > 800:
        windowHeight = 800

    windowWidth = int(round((1200/800) * windowHeight, -1))

    widget.resize(windowWidth, windowHeight)

    widget.show()
    sys.exit(app.exec_())
except Exception as e:
    logging.getLogger().info(str(e) + " main")
