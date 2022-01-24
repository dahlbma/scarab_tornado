# Pyinstaller how-to
Install Python 3.8 and pip, along with the requirements listed in 'requirements.txt'.

To compile main.py simply invoke
  
  `pyinstaller main.spec`

Make sure that you are using Python 3.8, along with the requirements listed in requirements.txt.
To make sure you are compiling with Python 3.8, try:

  `py -3.8 -m PyInstaller main.spec`

Which will create two directories: `build` \& `dist`, with `dist` containing the executable `chemreg`.

To install the launcher, use the same commands with 'main.spec' changed to 'launcher.spec'.
