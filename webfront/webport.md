# ChemReg Web Port Plan

## Recommendation

Use React with TypeScript and Vite for the browser UI. Keep Python, Tornado, RDKit, and the database as the system of record. The first phase should be a new client in `webfront/` that calls the existing HTTP endpoints; it should not rewrite chemistry logic, database access, or the schema.

## Current Baseline

- The desktop client is PyQt5. Its main workflows are sign-in, search and record browsing, compound registration/editing, bulk import, and metadata maintenance. See [the login screen](../chemReg/loginscreen.py), [search screen](../chemReg/searchscreen.py), [registration screen](../chemReg/regscreen.py), [bulk loader](../chemReg/sdfreg.py), and [metadata dialog](../chemReg/addmetatags.py).
- The existing HTTP client is [chemReg/dbInterface.py](../chemReg/dbInterface.py); Tornado registers its endpoints in [angular.py](../angular.py), with handlers in [dbInterface.py](../dbInterface.py).
- Useful existing endpoints include `/getDatabase`, `/login`, `/api/search`, `/api/getTextColumn`, `/api/getColComboData`, `/api/getMolfile`, `/api/loadMolfile`, `/api/createRegno`, `/api/update`, `/api/updateRegnoBatch`, `/api/bcpvsRegCompound`, and the bulk-import endpoints.
- Authentication has a legacy contract: `/login` returns JSON, and protected calls send that response body in a custom `token` header. Preserve this behind one frontend API adapter for phase 1; do not make UI components know about the format. The contract is implemented in [auth.py](../auth.py). Keep remote development and production traffic on the existing HTTPS `/chemreg/` endpoint; use plain HTTP only for a loopback Tornado instance.
- The Tornado server serves static files from `dist/`; the current [dist/index.html](../dist/index.html) is a launcher-download page. Do not overwrite it as part of local frontend development.
- Existing Qt forms use a blue-to-white sign-in background, a pale blue search screen (`#D4EFF8`), and pale red registration/import screens (`#FFCCCC`). Carry these recognizable color cues and dense form layout into responsive browser screens rather than copying the fixed desktop dimensions. References: [sign-in UI](../chemReg/assets/welcomescreen.ui), [search UI](../chemReg/assets/searchchem.ui), and [registration UI](../chemReg/assets/regchem.ui).

## Current Implementation

The initial React client is implemented in this directory and runs independently of the Python client. It includes:

- Sign-in using the legacy response-body `token` header contract, with session state kept in memory.
- Search by registration number, submitter, compound ID, and batch; result navigation, record details, molecular properties, and structure preview.
- Manual registration and editing using the existing draft, field-update, molfile-upload, and BCPVS endpoints.
- Supplier, library, and salt maintenance using the existing metadata endpoints.
- SDF, SMILES `.xlsx`/`.xls`, no-structure workbook, and CTRL compound import with field mapping, ELN allocation, progress, and a downloadable partial-results report.

No Tornado routes, Python dependencies, chemistry logic, or database schema have been changed. Before relying on write workflows, verify them against the Test database with authorized credentials; only parser, API-contract, build, and mocked UI checks have been run here.

## Phase 1: React Client, Existing Backend

1. **Confirm the first release scope.** Port sign-in, search/record details, manual registration/editing, bulk import, and metadata maintenance. Use the existing Test database for integration checks; do not test writes against Live.
2. **Create the frontend project in `webfront/`.** Set up React, TypeScript, and Vite. Keep its source and build output inside `webfront/`; do not change Python dependencies or backend routes.
3. **Add a single typed API client.** Match existing paths, query parameters, form-encoded requests, and multipart file uploads. Centralize response parsing for JSON, plain text, and non-200 responses. Configure Vite to proxy `/login`, `/getDatabase`, `/api/*`, and `/mols/*` through the existing HTTPS `/chemreg/` endpoint so browser requests are same-origin during development and no CORS change is needed.
4. **Implement sign-in and session state.** Load the Live/Test choices from `/getDatabase`, post credentials to `/login`, and keep the session in memory. The API client should preserve the exact login response text for the current `token` header contract, parse the response separately for display/state, and clear the session on logout or an authorization failure. Do not persist passwords or tokens in local storage.
5. **Build the app shell and visual tokens.** Match the existing screen colors, labels, field density, and navigation. Make layouts responsive, keyboard accessible, and explicit about loading, empty, validation, and error states.
6. **Port search and record browsing first.** Support the existing search fields (registration number, submitter, compound ID, and batch), hit count, previous/next navigation, record details, and molecule preview. Use existing endpoints and keep the initial slice read-oriented so it is easy to compare with the desktop client.
7. **Port manual registration and editing.** Reproduce the current metadata fields, combo-box data, batch validation, structure upload, molecule properties, and registration flow using existing endpoints. Preserve the current lifecycle and validation rules initially; avoid silently changing when edits are persisted.
8. **Port metadata maintenance.** Add the existing supplier, library, and salt workflows using their current endpoints and confirmation/validation behavior.
9. **Port bulk registration after the single-record flow works.** Support SDF, Excel SMILES, and no-structure Excel files, field mapping, ELN allocation, CTRL compounds, progress, and per-row failure reporting. Use existing per-record APIs initially and show clearly when an upload partially succeeds; the current backend does not provide an atomic batch job.
10. **Verify parity in the Test database.** Add frontend tests for the API adapter, form validation, search state, and import parsing. Smoke-test login, no-result and multi-result searches, navigation, valid/invalid edits, structure upload, each import type, ELN boundaries, CTRL handling, and partial failures. Compare the browser workflows with the existing Qt screens.
11. **Choose a staging deployment path without replacing the launcher page.** For phase 1, use Vite's proxy for development and a same-origin reverse-proxy path for a browser-accessible staging build. If the Tornado server must host the bundle directly, add a dedicated static path in a later deployment change; do not replace the current `dist/index.html` implicitly.

## Phase 2: Backend Follow-Up

After the React workflows are validated, make backend changes as a separate step:

1. Replace the legacy custom `token` header contract with a documented auth scheme and consistent session expiry/logout behavior. Add narrowly scoped CORS only if deployment requires a separate origin.
2. Standardize API responses and validation errors, then introduce versioned endpoints without breaking the desktop client during migration.
3. Add server-side validation and parameterized database operations for mutable/search inputs, with integration tests.
4. Consider a batch-import endpoint that validates before writing and returns a job/progress/result summary. This avoids many browser requests and makes partial failures easier to explain and recover from.
5. Finalize production static hosting and retain the existing launcher-download routes.

## Phase 1 Exit Criteria

- The web client uses the existing backend without schema or chemistry changes.
- Core desktop workflows have browser equivalents and pass against the Test database.
- Live writes remain disabled during development and test runs.
- The web build does not overwrite the launcher page or interfere with launcher downloads.
- Known compatibility limits (legacy auth header and non-atomic bulk import) are documented in the frontend API layer and user-facing error states.