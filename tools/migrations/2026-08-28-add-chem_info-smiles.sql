-- 2026-08-28: Store supplier SMILES separately from supplier molfiles.
-- Apply once on each environment via the mysql client:
--   mysql -h <host> -u <user> -p < 2026-08-28-add-chem_info-smiles.sql
--
-- This file is NOT idempotent. If the column already exists, comment out the
-- corresponding ALTER TABLE statement.

ALTER TABLE chem_reg.chem_info       ADD COLUMN SMILES mediumtext DEFAULT NULL;
ALTER TABLE chem_reg_test.chem_info  ADD COLUMN SMILES mediumtext DEFAULT NULL;
