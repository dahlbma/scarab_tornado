I have a compound registration program in the repo:https://github.com/dahlbma/scarab_tornado/
This repo is checked out in the directory: /home/mats.dahlberg/scarab_data/scarab_tornado

The file of interest in the repo is the dbInterface.py file.Read it and analyze the way structures and batches are registered.
Functions of special interest are: addStructure, registerNewCompound, registerNewBatch
Now to the problem. I have a file named chemreg_tofix_claude.csv. This is a CSV file with compounds that I would like to fix the structure for.It looks like:
compound_id	notebook_ref	source	regno	status	comment	mismatch_category	notes	registered_smiles	supplier_smiles
CBK004713	DT9037012	chem_reg	3938108	MISMATCH	CONNECTIVITY_DIFF	CONNECTIVITY_DIFF		NS(=O)(=O)c1ccc(NC(=O)c2cc3ccccc3o2)cc1	NS(=O)(=O)c1ccc(NC(=O)C2Cc3ccccc3O2)cc1
CBK005606	EE5502872	chem_reg	4004661	MISMATCH	CONNECTIVITY_DIFF	CONNECTIVITY_DIFF		O=C(Nc1cccnc1)C1CCCS1	O=C(Nc1cccnc1)c1cccs1

Basically we have a mysql table "chem_reg.chem_info" that contains the data we get from the supplier of the compounds. What I would like to do is to make sure that the structure in the column: chem_reg.chem_info.MOLFILE represents the same structure as in the bcpvs.compound table (smiles_std_smiles_std_string, inchi, inchi_key, mf, SEP_MOL_MONOISO_MASS).
But we have to be very careful, we can't simply update the structure elements in the bcpvs.compound table, because there might be other batches in the batch table: bcpvs.batch that points to the same structure and if we blindly updates the structure these other batches might get the wrong structure. So in the cases where there are multiple batches pointing to the structure we want to update we must go back and verify the structure of all the batches. These batches can have their original structure in either the chem_reg.chem_info table or in the chemspec.chem_info table. You can tell what table to look in by checking the bcpvs.batch.chemreg_regno or bcpvs.batch.chemspec_regno. Only one of these two columns will have a value and that tells if the batch was registered in the old chemspec system or in the new chem_reg system.
In the easy case there are no other batches pointing to the compound_id or all the batches pointing to the compound_id we want to update have the same structure in chem_reg and chemspec, then we can just update the bcpvs.compound table with the correct structure information. Also notice that we use the molsoft cartridge (molcart) for chemical searches, so any compound_ids that gets an updated structure in bcpvs.compound we also need to update the tables that molsoft uses, how this is done you can check in the function addStructure() in the dbInterface.py file in the scarab_tornado repo.
In the more complicated case where there are multiple batches and they have different structures in the chemspec and/or chem_reg systems we might have to register a new compound_id with the correct structure in the bcpvs.compound table.


We use a chemical cartridge named molcart. The tables that molcart uses to store structure information are: bcpvs.JCMOL_MOLTABLE, bcpvs.JCMOL_MOLTABLE_ukey, bcpvs.JCMOL_MOLTABLE_MOL_keysim, bcpvs.JCMOL_MOLTABLE_MOL_key

I have created a python file 'update_structures.py' that reads the chemreg_tofix_claude.csv that investigates the problem. The aim is to update the structures in bcpvs.compound and the molcart tables But for now I would like you to run in "dry mode", DO NOT update anything in the database yet! We need to be very careful and first generate a file that tells me what updates are going to be done. This is so that I can read and go through the proposed changes and make sure they are safe and correct.

We use the python RDKit library and the version installed on the server is 2022.3. We have created a set of python scripts that tries to pinpoint what the cause of the error is. I first thought that it was the RDKit library that was the cause of the problem. But the test scripts have not been able to prove this. These test scripts are:
test_smiles_provenance.py
test_server_molfile.py
test_server_env.py
test_sdf_compare.py
test_mysql_roundtrip.py
test_final_diagnostic.py
test_cross_process.py
test_tautomer_analysis.py

None of these tests has been able to pinpoint where the error is introduced.

I registered 2 sdfiles (CBCS0523_CC04640_192cpds_Carlsson.sdf and CBCS0523_CC04606_208cpds_Carlsson.sdf) that contains the same structures and when I registered the second file there was 11 structures that got registered as new compounds even though they where present in the first file. Here are these compounds:
compound_id	notebook_ref	source	regno	status	comment	mismatch_category	notes	registered_smiles	supplier_smiles
CBK028998	EH8331120	chem_reg	4052394	MISMATCH	CONNECTIVITY_DIFF	CONNECTIVITY_DIFF		O=C(NCC1CCCO1)NC1CCCCC1	O=C(NCC1CCCO1)Nc1ccccc1
CBK228517	EH8331126	chem_reg	4052400	MISMATCH	CONNECTIVITY_DIFF	CONNECTIVITY_DIFF		O=C(Nc1ccccc1)c1ccccc1	O=C(NC1CCCCC1)c1ccccc1
CBK278468	DT9019206	chem_reg	3919587	MISMATCH	CONNECTIVITY_DIFF	CONNECTIVITY_DIFF		NC(=O)c1ccncc1	NC(=O)C1CCNCC1
CBK278468	EH8331174	chem_reg	4052448	MISMATCH	CONNECTIVITY_DIFF	CONNECTIVITY_DIFF		NC(=O)c1ccncc1	NC(=O)C1CCNCC1
CBK278635	EH8331011	chem_reg	4052285	MISMATCH	CONNECTIVITY_DIFF	CONNECTIVITY_DIFF		O=C(Nc1ccccc1)c1cccc(O)c1	O=C(NC1CCCCC1)c1cccc(O)c1
CBK291635	EH8331059	chem_reg	4052333	MISMATCH	CONNECTIVITY_DIFF	CONNECTIVITY_DIFF		O=C(Nc1ccccc1)c1cccnc1	O=C(NC1CCCCC1)c1cccnc1
CBK291792	EH8331113	chem_reg	4052387	MISMATCH	CONNECTIVITY_DIFF	CONNECTIVITY_DIFF		CNC(=O)C1CCNCC1	CNC(=O)c1ccncc1
CBK389719	EH8331149	chem_reg	4052423	MISMATCH	CONNECTIVITY_DIFF	CONNECTIVITY_DIFF		O=C(NC1CCNCC1)C1CCCCC1	O=C(NC1CCNCC1)c1ccccc1
CBK699880	EH8331026	chem_reg	4052300	MISMATCH	CONNECTIVITY_DIFF	CONNECTIVITY_DIFF		O=C(NCc1ccccc1)c1ccccc1	O=C(NCc1ccccc1)C1CCCCC1
CBK714175	EH8331106	chem_reg	4052380	MISMATCH	CONNECTIVITY_DIFF	CONNECTIVITY_DIFF		CNC(=O)c1cccnc1	CNC(=O)C1CCCNC1
CBK714200	EH8331132	chem_reg	4052406	MISMATCH	CONNECTIVITY_DIFF	CONNECTIVITY_DIFF		O=C(NCC1CCCO1)c1cccc(Cl)c1	O=C(NCc1ccco1)c1cccc(Cl)c1


I want to introduce logging information in db.Interface.py that catches information when a structure is about to get registered in the case that there is a chemical difference between the structure entered into the chem_reg.chem_info.MOLFILE and the structure of the bcpvs.compound_smiles_std.
Also make sure that when there is a difference in the structure between chem_reg.chem_info.MOLFILE and the SMILES that is about to be entered I want to procedure to not insert the new compound into bcpvs.compound. The user that are uploading the sdfile through must be notified that the structure failed. All compounds that fail to register during a sdfile registrtion (through the sdfreg.py) needs to be stored in a sdfile on the users computer so the user can examine the problematic entries.




The tables look like:
CREATE TABLE bcpvs.batch (
  `compound_id` varchar(26) CHARACTER SET utf8mb3 COLLATE utf8mb3_bin NOT NULL,
  `notebook_ref` varchar(16) NOT NULL,
  `submitter` varchar(30) DEFAULT NULL,
  `SUBMITTAL_DATE` datetime DEFAULT NULL,
  `project` varchar(30) DEFAULT NULL,
  `PURITY` double DEFAULT NULL,
  `supplier` varchar(100) DEFAULT NULL,
  `AMOUNT` double DEFAULT NULL,
  `BATCH_COMMENT` text,
  `BIOLOGICAL_MW` double DEFAULT NULL,
  `colour` varchar(20) DEFAULT NULL,
  `physical_state` varchar(20) DEFAULT NULL,
  `MELT_PT_LOW` double DEFAULT NULL,
  `MELT_PT_HIGH` text,
  `cryst_solvent` varchar(20) DEFAULT NULL,
  `CRYST_SOLVENT_RATIO` double DEFAULT NULL,
  `CREATED_DATE` datetime DEFAULT NULL,
  `CREATED_BY_PKEY` text,
  `UPDATED_DATE` text,
  `UPDATED_BY_PKEY` text,
  `library_id` varchar(100) DEFAULT NULL,
  `compound_type` varchar(20) DEFAULT NULL,
  `product_type` varchar(20) DEFAULT NULL,
  `supplier_id` varchar(200) DEFAULT NULL,
  `BEE_PURITY` double DEFAULT NULL,
  `MEASURED_MONOISOMASS` double DEFAULT NULL,
  `supplier_batch` varchar(100) DEFAULT NULL,
  `chemspec_regno` bigint DEFAULT NULL,
  `suffix` text,
  `chemreg_regno` int DEFAULT NULL,
  `restriction_comment` varchar(45) DEFAULT NULL,
  `oldcbk` varchar(26) DEFAULT NULL,
  PRIMARY KEY (`notebook_ref`),
  UNIQUE KEY `batch_idx` (`notebook_ref`),
  KEY `compound_id` (`compound_id`),
  KEY `chems_regn_idx` (`chemspec_regno`),
  KEY `chemreg_regno` (`chemreg_regno`),
  KEY `fk_batch_library` (`library_id`),
  CONSTRAINT `batch_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `compound` (`compound_id`),
  CONSTRAINT `batch_ibfk_2` FOREIGN KEY (`chemreg_regno`) REFERENCES `chem_reg`.`chem_info` (`regno`),
  CONSTRAINT `batch_ibfk_3` FOREIGN KEY (`chemspec_regno`) REFERENCES `chemspec`.`chem_info` (`regno`),
  CONSTRAINT `fk_batch_library` FOREIGN KEY (`library_id`) REFERENCES `compound_library` (`library_name`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

CREATE TABLE bcpvs.compound (
  `compound_id` varchar(26) CHARACTER SET utf8mb3 COLLATE utf8mb3_bin DEFAULT NULL,
  `COMPOUND_ID_NUMERIC` bigint DEFAULT NULL,
  `suffix` text,
  `STEREO_COMMENTS` text,
  `COMPOUND_NAME` text,
  `STRUCTURE_COMMENTS` text,
  `CAS_NUMBER` text,
  `BEILSTEIN_REG_NO` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `CREATED_BY_PKEY` text,
  `UPDATED_DATE` text,
  `UPDATED_BY_PKEY` text,
  `mf` varchar(200) DEFAULT NULL,
  `CAS_SEARCH_DATE` datetime DEFAULT NULL,
  `BEILSTEIN_SEARCH_DATE` datetime DEFAULT NULL,
  `ip_rights` varchar(30) DEFAULT NULL,
  `SEP_MOL_MONOISO_MASS` double DEFAULT NULL,
  `smiles_std` varchar(1000) CHARACTER SET utf8mb4 COLLATE utf8mb4_bin DEFAULT NULL,
  `smiles_std_notaut` varchar(1000) DEFAULT NULL,
  `smiles_std_string` varchar(1000) CHARACTER SET utf8mb4 COLLATE utf8mb4_bin DEFAULT NULL,
  `chembl_id` varchar(20) CHARACTER SET utf8mb3 COLLATE utf8mb3_general_ci DEFAULT NULL,
  `inchi` varchar(4000) DEFAULT NULL,
  `inchi_key` varchar(27) DEFAULT NULL,
  `pref_name` varchar(255) DEFAULT NULL,
  UNIQUE KEY `compound_idx` (`compound_id`),
  KEY `idnumeric_idx` (`COMPOUND_ID_NUMERIC`),
  KEY `idx_smiles_std` (`smiles_std`(255)),
  KEY `idx_inchi_key` (`inchi_key`),
  KEY `idx_chembl_id` (`chembl_id`),
  CONSTRAINT `fk_compound_chembl_id` FOREIGN KEY (`chembl_id`) REFERENCES `chembl_35`.`molecule_dictionary` (`chembl_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;


CREATE TABLE chem_reg.chem_info (
  `regno` int DEFAULT NULL,
  `chemist` varchar(40) DEFAULT NULL,
  `JPAGE` varchar(16) DEFAULT NULL,
  `COMPOUND_ID` varchar(16) DEFAULT NULL,
  `RDATE` datetime DEFAULT NULL,
  `COMMENTS` mediumtext,
  `EXTERNAL_ID` varchar(200) DEFAULT NULL,
  `project` varchar(40) DEFAULT NULL,
  `PRODUCT` varchar(30) DEFAULT NULL,
  `LIBRARY_ID` varchar(100) DEFAULT NULL,
  `SOURCE` varchar(200) DEFAULT NULL,
  `COMPOUND_TYPE` varchar(20) DEFAULT NULL,
  `REG_DATE` datetime DEFAULT NULL,
  `CHEM_NAME` mediumtext,
  `CREATED_BY_PKEY` double DEFAULT NULL,
  `SUPPLIER_BATCH` mediumtext,
  `C_MF` mediumtext,
  `C_MW` double DEFAULT NULL,
  `C_MONOISO` double DEFAULT NULL,
  `C_CHNS` mediumtext,
  `MOLFILE` mediumtext,
  `CHROM_TEXT` mediumtext,
  `NMR_TEXT` mediumtext,
  `MS_TEXT` mediumtext,
  `CHROM_PURITY` int DEFAULT NULL,
  `MS_PURITY` int DEFAULT NULL,
  `NMR_PURITY` int DEFAULT NULL,
  `SOLVENT` varchar(40) DEFAULT NULL,
  `PURITY` int DEFAULT NULL,
  `suffix` text,
  `IP_RIGHTS` varchar(20) DEFAULT NULL,
  `sdfile_sequence` int DEFAULT NULL,
  UNIQUE KEY `regno_pk` (`regno`),
  UNIQUE KEY `index_name` (`JPAGE`),
  KEY `cmp_idx` (`COMPOUND_ID`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;


CREATE TABLE chemspec.chem_info (
  `regno` bigint DEFAULT NULL,
  `CHEMIST` text,
  `jpage` varchar(20) DEFAULT NULL,
  `RDATE` datetime DEFAULT NULL,
  `QUANTITY` text,
  `solvent` varchar(50) DEFAULT NULL,
  `STORAGE` text,
  `COMMENTS` text,
  `EXTERNAL_ID` text,
  `project` varchar(40) DEFAULT NULL,
  `ACCMS_ANALYSIS` text,
  `RETURN_SAMPLE` text,
  `STR_COMMENTS` text,
  `HNMR_ANALYSIS` text,
  `CNMR_ANALYSIS` text,
  `MS_ANALYSIS` text,
  `CHNS_ANALYSIS` text,
  `ADV_ANALYSIS` text,
  `PRODUCT` text,
  `CHNS_COMMENT` text,
  `PURITY` double DEFAULT NULL,
  `ID_OK` text,
  `INTERP_REQUEST` text,
  `CHROM_ANALYSIS` text,
  `library_id` varchar(80) DEFAULT NULL,
  `SOURCE` text,
  `PLATE_ID` text,
  `CONFIG_ID` text,
  `COMPOUND_TYPE` text,
  `REG_DATE` datetime DEFAULT NULL,
  `CHEM_NAME` text,
  `OK_FOR_REG` text,
  `OK_FOR_REG_SIGN` text,
  `STEREO_COMMENTS` text,
  `CREATED_BY_PKEY` double DEFAULT NULL,
  `CREATED_DATE` text,
  `UPDATED_BY_PKEY` text,
  `UPDATED_DATE` text,
  `IR_ANALYSIS` text,
  `OTHER_ANALYSIS` text,
  `SUPPLIER_BATCH` text,
  `C_MF` text,
  `C_MW` double DEFAULT NULL,
  `C_MONOISO` double DEFAULT NULL,
  `C_CHNS` text,
  `EXPERIMENTAL` text,
  `COLOUR` text,
  `PHYSICAL_STATE` text,
  `molfile` text,
  UNIQUE KEY `regno_pk` (`regno`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
