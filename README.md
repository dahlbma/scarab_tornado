Backend for chemical registration application is using rdkit for chemical functions.
Frontend in python using pyqt5.


The backend needs a config file config.py with content:

```
database = dict(
    host = 'hostname',
    user = 'username',
    password = 'password'
)

secret_key = 'Some random string'
```
