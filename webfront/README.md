# ChemReg Web Frontend

The React client is kept separate from Tornado. During development, Vite serves the UI on port `8082` and proxies API and molecule-image requests to the existing HTTPS ChemReg endpoint at `https://esox3.scilifelab.se/chemreg/`.

Current screens cover sign-in, compound search and details, manual registration/editing, supplier/library/salt maintenance, and bulk registration. Bulk import accepts SDF files, `.xlsx`/`.xls` workbooks with a `Smile` column, and workbooks with the exact no-structure columns from the desktop workflow. The app keeps the existing Tornado API and database as the system of record.

```sh
npm install
npm run dev
```

Run the frontend checks with:

```sh
npm test
npm run build
npm audit
```

To use a local Tornado instance or a different backend address, set `CHEMREG_API_TARGET` when starting Vite:

```sh
CHEMREG_API_TARGET=http://127.0.0.1:8081 npm run dev
```

Keep production or remote targets on HTTPS because sign-in credentials and the legacy token header traverse this proxy.

Use the Test database while verifying write workflows. Bulk registration posts rows individually; if a later row fails, earlier successful rows remain registered and are listed in the import report.