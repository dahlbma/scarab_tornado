import sys, os, logging, traceback
from PyQt5.uic import loadUi
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QApplication, QDialog, QMainWindow
from PyQt5 import QtGui

from chemreglib import *

def error_handler(etype, value, tb):
    err_msg = "".join(traceback.format_exception(etype, value, tb))
    logger.exception(err_msg)

class LauncherScreen(QDialog):
    def __init__(self):
        super(LauncherScreen, self).__init__()
        self.mod_name = "launcher"
        logger = logging.getLogger(self.mod_name)
        loadUi(resource_path("assets/launcher.ui"), self)
        self.update_chemreg_btn.clicked.connect(self.updatefunction)
        self.update_chemreg_btn.setDefault(False)
        self.update_chemreg_btn.setAutoDefault(False)
        self.run_chemreg_btn.clicked.connect(self.runfunction)
        self.run_chemreg_btn.setDefault(True)
        self.run_chemreg_btn.setAutoDefault(True)

        if self.ver_check() is False: # outdated
            self.status_lab.setText("""ChemReg is outdated!<br>
            Please <b>'Update ChemReg'</b> or<br>
            <b>'Run ChemReg'</b> to Update.""")


    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Return or event.key() == QtCore.Qt.Key_Enter:
            self.runfunction()
        elif event.key() == QtCore.Qt.Key_R: 
            self.updatefunction()

    def ver_check(self):
        # return true if chemReg is outdated
        try: 
            r = requests.get('<>') # get file version
        except Exception as e:
            self.status_lab.setText("ERROR no connection")
            logging.getLogger(self.mod_name).error(str(e))
            return
        info_dict = dict()
        with open('./ver.dat') as f:
            for line in f:
                # reads each line and trims of extra the spaces 
                # and gives only the valid words
                command, description = line.strip().split(None, 1)
                info_dict[command] = description.strip()
        # check if versions match
        ok = True if r.content == info_dict['version'] else False
        return ok


    def updatefunction(self):
        # check if versions match
        match = self.ver_check()
        
        if match:
            # all is well
            return 0
        else:
            # update
            os_name = platform.system()
            try: 
                r = requests.get('<>', data={'os_name':os_name}) # fetch chemreg dist
            except Exception as e:
                self.status_lab.setText("ERROR ")
                logging.getLogger(self.mod_name).error(str(e))
                return


    def runfunction(self):
        self.updatefunction()
        os_name = platform.system()
        if os_name == 'Windows':
            subprocess.Popen(['chemreg.exe'], shell=True)
        elif os_name == 'Linux':
            subprocess.Popen(['./chemreg'], shell=True)
        elif os_name == 'Darwin':
            subprocess.Popen(['open', 'chemreg'], shell=True)
        else:
            send_msg("Error", "Can not launch chemreg, unknown OS", icon=QMessageBox.Warning)
        QtWidgets.QApplication.instance().quit()
        return
        

#base settings for logging
level=logging.INFO

# init root logger
logger = logging.getLogger()
logger.setLevel(level)

# console logging
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter('%(name)s:%(message)s')
ch.setFormatter(formatter)

# file logging
file=os.path.join(".","chemreg_launcher.log")
fh = logging.FileHandler(file)
fh.setLevel(level)
formatter = logging.Formatter('%(asctime)s : %(name)s:%(levelname)s : %(message)s',
                              datefmt='%m/%d/%Y %H:%M:%S')
fh.setFormatter(formatter)

logger.addHandler(ch)
logger.addHandler(fh)


try:
    # base app settings
    os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "2"
    app = QApplication(['Chem Reg Launcher'])
    clipboard = app.clipboard()

    launch = LauncherScreen()
    widget = QtWidgets.QStackedWidget()
    widget.addWidget(launch)

    desktop = QApplication.desktop()
    windowHeight = 340
    windowWidth = 508

    #windowHeight = int(round(0.5 * desktop.screenGeometry().height(), -1))
    #if windowHeight > 800:
    #    windowHeight = 800

    #windowWidth = int(round((1200/800) * windowHeight, -1))

    widget.resize(windowWidth, windowHeight)

    widget.show()
    app.setWindowIcon(QtGui.QIcon('asssets/chem.ico'))
    widget.setWindowIcon(QtGui.QIcon('assets/chem.ico'))
    sys.exit(app.exec_())
except Exception as e:
    logger.info(str(e))
