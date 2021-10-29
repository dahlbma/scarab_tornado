# Pyinstaller how-to
Install python, along with the following packages using pip or your favourite package manager:
* PyQt5
* requests
* pyinstaller

Now, to compile main.py simply invoke
  
  `pyinstaller main.spec`

Which will create two directories: `build` \& `dist`, with `dist` containing the executable `chemreg`.
