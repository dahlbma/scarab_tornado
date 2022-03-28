Backend for chemical registration application is using rdkit for chemical functions.
Frontend in python using pyqt5.

## Frontend
### PyInstaller How-To
Currently only builds on `Python 3.8`, with required package versions listed in `requirements.txt`.
With `frontend` as current directory, build the main ChemReg executable with:

    <pyinstaller> main.spec
which will build the main executable `ch`(.exe)
</br>or

    <pyinstaller> launcher.spec
which will build the launcher executable `chemreg`(.exe).

Substitute `<pyinstaller>` with your local appropriate PyInstaller module (possibly `py -3.8 -m PyInstaller` or just `python3 pyinstaller`, case sensitive module names).

<!--
<b> Make sure to build a new version of the `ch` executable before using `upload.py` to upload a new version!
</b>
-->

## Backend
The backend needs a config file config.py with content:

```
database = dict(
    host = 'hostname',
    user = 'username',
    password = 'password'
)

secret_key = 'Some random string'
```