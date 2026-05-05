-- 2026-04-29: Propagate the standardized InChIKey from ChemRegAddMol into
-- chem_info so BcpvsRegCompound does not have to re-run RDKit's
-- TautomerEnumerator.Canonicalize() on the supplier molfile (which is
-- non-deterministic in RDKit 2022.3 for some scaffolds and is the root
-- cause of the duplicate-compound rows the structureCorrection tools
-- have to repair).
--
-- Apply once on each environment via the mysql client:
--   mysql -h <host> -u <user> -p < 2026-04-29-add-chem_info-inchi_key.sql
--
-- This file is NOT idempotent. If the column or index already exists,
-- comment out the corresponding line. (For an idempotent backfill use
-- backfill_inchi_key.py.)

ALTER TABLE chem_reg.chem_info       ADD COLUMN inchi_key VARCHAR(27) DEFAULT NULL;
CREATE INDEX idx_chem_info_inchi_key ON chem_reg.chem_info (inchi_key);

ALTER TABLE chem_reg_test.chem_info  ADD COLUMN inchi_key VARCHAR(27) DEFAULT NULL;
CREATE INDEX idx_chem_info_inchi_key ON chem_reg_test.chem_info (inchi_key);
