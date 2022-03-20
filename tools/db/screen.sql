-- MySQL dump 10.13  Distrib 8.0.26, for Linux (x86_64)
--
-- Host: localhost    Database: screen
-- ------------------------------------------------------
-- Server version	8.0.26

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!50503 SET NAMES utf8mb4 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `aa_cyp3a4_stab`
--

DROP TABLE IF EXISTS `aa_cyp3a4_stab`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `aa_cyp3a4_stab` (
  `PKEY` bigint DEFAULT NULL,
  `ABASE_TTORDIRNO` bigint DEFAULT NULL,
  `ABASE_TTOROBJSEQNO` bigint DEFAULT NULL,
  `project` varchar(20) DEFAULT NULL,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `compound_batch` varchar(20) DEFAULT NULL,
  `TDATE` datetime DEFAULT NULL,
  `STABILITY` text,
  `ERROR_FLAG` text,
  `PCTSTART` text,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `ABASE_TOCCID` text,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `aa_cyp3a4_stab_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `aa_elph_1`
--

DROP TABLE IF EXISTS `aa_elph_1`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `aa_elph_1` (
  `PKEY` bigint DEFAULT NULL,
  `ABASE_TTORDIRNO` double DEFAULT NULL,
  `ABASE_TTOROBJSEQNO` double DEFAULT NULL,
  `project` varchar(20) DEFAULT NULL,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `compound_batch` varchar(20) DEFAULT NULL,
  `PROBE_METHOD` text,
  `TDATE` datetime DEFAULT NULL,
  `RESULT` text,
  `RELFLUOR` double DEFAULT NULL,
  `ERROR_FLAG` text,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `MDATE` datetime DEFAULT NULL,
  `ABASE_TOCCID` text,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `aa_elph_1_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `aa_ez4u_screen`
--

DROP TABLE IF EXISTS `aa_ez4u_screen`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `aa_ez4u_screen` (
  `PKEY` bigint DEFAULT NULL,
  `ABASE_TTORDIRNO` bigint DEFAULT NULL,
  `ABASE_TTOROBJSEQNO` bigint DEFAULT NULL,
  `PROJECT` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `COMPOUND_BATCH` text,
  `TDATE` datetime DEFAULT NULL,
  `AGE` double DEFAULT NULL,
  `RESULT` text,
  `TC50` double DEFAULT NULL,
  `EXPOSURE_TIME` double DEFAULT NULL,
  `CELL_LINE` text,
  `ERROR_FLAG` text,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `MDATE` datetime DEFAULT NULL,
  `CONC` text,
  `ABASE_TOCCID` text,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `aa_ez4u_screen_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `aa_induction`
--

DROP TABLE IF EXISTS `aa_induction`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `aa_induction` (
  `PKEY` bigint DEFAULT NULL,
  `ABASE_TTORDIRNO` bigint DEFAULT NULL,
  `ABASE_TTOROBJSEQNO` bigint DEFAULT NULL,
  `ABASE_TOCCID` text,
  `PROJECT` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `COMPOUND_BATCH` text,
  `TDATE` datetime DEFAULT NULL,
  `ENZYME` text,
  `CELL_LINE` text,
  `METHOD` text,
  `CONC` double DEFAULT NULL,
  `INDUCTION` double DEFAULT NULL,
  `EC50` double DEFAULT NULL,
  `REF_COMPOUND` text,
  `REF_INDUCTION` double DEFAULT NULL,
  `RATIO` double DEFAULT NULL,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `MDATE` datetime DEFAULT NULL,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `aa_induction_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `aa_ips`
--

DROP TABLE IF EXISTS `aa_ips`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `aa_ips` (
  `PKEY` bigint DEFAULT NULL,
  `ABASE_TTORDIRNO` double DEFAULT NULL,
  `ABASE_TTOROBJSEQNO` double DEFAULT NULL,
  `PROJECT` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `compound_batch` varchar(16) DEFAULT NULL,
  `TDATE` datetime DEFAULT NULL,
  `AREAPCT` double DEFAULT NULL,
  `ION_MODE` text,
  `ISOMERS` text,
  `PURITY_NUMERIC` double DEFAULT NULL,
  `SPURITY` double DEFAULT NULL,
  `STABILITY` double DEFAULT NULL,
  `IDENTITY` text,
  `ELUTES_IN_FRONT` text,
  `NO_MS` text,
  `NO_UV` text,
  `PURITY` text,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `MDATE` datetime DEFAULT NULL,
  `RAW_DATA_FILE` text,
  `ABASE_TOCCID` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `RECORD_SOURCE` text,
  `RECORD_SOURCE_PKEY` text,
  KEY `compound_batch` (`compound_batch`),
  CONSTRAINT `aa_ips_ibfk_1` FOREIGN KEY (`compound_batch`) REFERENCES `bcpvs`.`batch` (`notebook_ref`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `aa_ltssolubility`
--

DROP TABLE IF EXISTS `aa_ltssolubility`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `aa_ltssolubility` (
  `PKEY` bigint DEFAULT NULL,
  `ABASE_TTORDIRNO` bigint DEFAULT NULL,
  `ABASE_TTOROBJSEQNO` bigint DEFAULT NULL,
  `ABASE_TOCCID` text,
  `PROJECT` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `COMPOUND_BATCH` text,
  `TDATE` datetime DEFAULT NULL,
  `METHOD` text,
  `BUFFER` text,
  `IONIC_STRENGTH` double DEFAULT NULL,
  `SOL_MASS` double DEFAULT NULL,
  `SOL_MOL` double DEFAULT NULL,
  `TEMP` text,
  `TIME` bigint DEFAULT NULL,
  `PH` double DEFAULT NULL,
  `PH2` double DEFAULT NULL,
  `R_SQ` double DEFAULT NULL,
  `DATA_QUALITY` text,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `MDATE` datetime DEFAULT NULL,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `aa_ltssolubility_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `aa_metstab`
--

DROP TABLE IF EXISTS `aa_metstab`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `aa_metstab` (
  `PKEY` bigint DEFAULT NULL,
  `ABASE_TTORDIRNO` double DEFAULT NULL,
  `ABASE_TTOROBJSEQNO` double DEFAULT NULL,
  `PROJECT` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `COMPOUND_BATCH` text,
  `TDATE` datetime DEFAULT NULL,
  `GENDER` text,
  `CLH` double DEFAULT NULL,
  `CLINT_TEST_TUBE` double DEFAULT NULL,
  `CLINT` double DEFAULT NULL,
  `CLINT2` text,
  `T_HALF` double DEFAULT NULL,
  `E` double DEFAULT NULL,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `METHOD` text,
  `SPECIES` text,
  `STRAIN` text,
  `MDATE` datetime DEFAULT NULL,
  `ABASE_TOCCID` text,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `aa_metstab_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `aa_mmt_ames`
--

DROP TABLE IF EXISTS `aa_mmt_ames`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `aa_mmt_ames` (
  `PKEY` bigint DEFAULT NULL,
  `ABASE_TTORDIRNO` bigint DEFAULT NULL,
  `ABASE_TTOROBJSEQNO` bigint DEFAULT NULL,
  `ABASE_TOCCID` text,
  `PROJECT` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `COMPOUND_BATCH` text,
  `TDATE` datetime DEFAULT NULL,
  `STRAIN` text,
  `METHOD` text,
  `MUTAGENICITY` text,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `MDATE` datetime DEFAULT NULL,
  `DATA_QUALITY` text,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `aa_mmt_ames_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `aa_p450inhib`
--

DROP TABLE IF EXISTS `aa_p450inhib`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `aa_p450inhib` (
  `PKEY` bigint DEFAULT NULL,
  `ABASE_TTORDIRNO` bigint DEFAULT NULL,
  `ABASE_TTOROBJSEQNO` bigint DEFAULT NULL,
  `ABASE_TOCCID` text,
  `PROJECT` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `COMPOUND_BATCH` text,
  `TDATE` datetime DEFAULT NULL,
  `ENZYME` text,
  `METHOD` text,
  `CONC` double DEFAULT NULL,
  `INHIBITION` double DEFAULT NULL,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `MDATE` datetime DEFAULT NULL,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `aa_p450inhib_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `aa_permeabilitet`
--

DROP TABLE IF EXISTS `aa_permeabilitet`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `aa_permeabilitet` (
  `PKEY` bigint DEFAULT NULL,
  `ABASE_TTORDIRNO` double DEFAULT NULL,
  `ABASE_TTOROBJSEQNO` double DEFAULT NULL,
  `PROJECT` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `COMPOUND_BATCH` text,
  `TDATE` datetime DEFAULT NULL,
  `CONC` double DEFAULT NULL,
  `EFFLUX_RATIO` double DEFAULT NULL,
  `MASS_BAL_AB` double DEFAULT NULL,
  `MASS_BAL_BA` double DEFAULT NULL,
  `PERM_CLASS` text,
  `ERROR_FLAG` text,
  `PAPP_AB` double DEFAULT NULL,
  `PAPP_BA` double DEFAULT NULL,
  `REF_CMPD1_ID` text,
  `REF_CMPD1_PAPP` double DEFAULT NULL,
  `REF_CMPD2_ID` text,
  `REF_CMPD2_PAPP` double DEFAULT NULL,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `MDATE` datetime DEFAULT NULL,
  `ABASE_TOCCID` text,
  `ASSAY_TYPE` text,
  `CELL_LINE` text,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `aa_permeabilitet_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `aa_pka`
--

DROP TABLE IF EXISTS `aa_pka`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `aa_pka` (
  `PKEY` bigint DEFAULT NULL,
  `ABASE_TTORDIRNO` bigint DEFAULT NULL,
  `ABASE_TTOROBJSEQNO` bigint DEFAULT NULL,
  `PROJECT` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `COMPOUND_BATCH` text,
  `TDATE` datetime DEFAULT NULL,
  `IONIC_STRENGTH` double DEFAULT NULL,
  `LOG_KD` double DEFAULT NULL,
  `METHOD` text,
  `TEMPERATURE` bigint DEFAULT NULL,
  `PH_RANGE` text,
  `PKA1_ACIDIC` double DEFAULT NULL,
  `PKA2_ACIDIC` double DEFAULT NULL,
  `PKA3_ACIDIC` text,
  `PKA1_ALKALINE` double DEFAULT NULL,
  `PKA2_ALKALINE` double DEFAULT NULL,
  `PKA3_ALKALINE` double DEFAULT NULL,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `MDATE` datetime DEFAULT NULL,
  `ABASE_TOCCID` text,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `aa_pka_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `aa_protein_binding`
--

DROP TABLE IF EXISTS `aa_protein_binding`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `aa_protein_binding` (
  `PKEY` bigint DEFAULT NULL,
  `ABASE_TTORDIRNO` bigint DEFAULT NULL,
  `ABASE_TTOROBJSEQNO` bigint DEFAULT NULL,
  `PROJECT` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `COMPOUND_BATCH` text,
  `TDATE` datetime DEFAULT NULL,
  `METHOD` text,
  `BINDING` text,
  `CONC` double DEFAULT NULL,
  `STABILITY` double DEFAULT NULL,
  `STABILITY24H` double DEFAULT NULL,
  `FB` double DEFAULT NULL,
  `FU` double DEFAULT NULL,
  `ERROR_FLAG` text,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `SPECIES` text,
  `N` text,
  `STABLE` text,
  `MDATE` datetime DEFAULT NULL,
  `ABASE_TOCCID` text,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `aa_protein_binding_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `aa_solubility`
--

DROP TABLE IF EXISTS `aa_solubility`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `aa_solubility` (
  `PKEY` bigint DEFAULT NULL,
  `ABASE_TTORDIRNO` double DEFAULT NULL,
  `ABASE_TTOROBJSEQNO` double DEFAULT NULL,
  `PROJECT` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `compound_batch` varchar(16) DEFAULT NULL,
  `TDATE` datetime DEFAULT NULL,
  `SOL_MASS` text,
  `BUFFER` text,
  `MEDIA_ADDITIVES` text,
  `SOL_MOL` double DEFAULT NULL,
  `ERROR_FLAG` text,
  `STOCK_CONC` text,
  `STOCK_SOLVENT` text,
  `TEMPERATURE` text,
  `PH` text,
  `R_SQUARED` text,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `MDATE` datetime DEFAULT NULL,
  `ABASE_TOCCID` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `RECORD_SOURCE` text,
  `RECORD_SOURCE_PKEY` text,
  KEY `compound_batch` (`compound_batch`),
  CONSTRAINT `aa_solubility_ibfk_1` FOREIGN KEY (`compound_batch`) REFERENCES `bcpvs`.`batch` (`notebook_ref`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `aa_spec_analysis`
--

DROP TABLE IF EXISTS `aa_spec_analysis`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `aa_spec_analysis` (
  `PKEY` text,
  `ABASE_TTORDIRNO` text,
  `ABASE_TTOROBJSEQNO` text,
  `PROJECT` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `COMPOUND_BATCH` text,
  `TDATE` text,
  `AREAPCT` text,
  `CHIRAL_PURITY` text,
  `FORMULATION` text,
  `IDENTITY` text,
  `ION_MODE` text,
  `ISOMERS` text,
  `PURITY` text,
  `SPURITY` text,
  `STABILITY` text,
  `SAMPLE_CONC` text,
  `SAMPLE_CONC2` text,
  `TEMPERATURE` text,
  `TIME` text,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `MDATE` text,
  `ABASE_TOCCID` text,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `aa_spec_analysis_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ba_binding_sp`
--

DROP TABLE IF EXISTS `ba_binding_sp`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ba_binding_sp` (
  `PKEY` bigint DEFAULT NULL,
  `ABASE_TTORDIRNO` bigint DEFAULT NULL,
  `ABASE_TTOROBJSEQNO` bigint DEFAULT NULL,
  `ABASE_TOCCID` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `PROJECT` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `COMPOUND_BATCH` text,
  `METHOD` text,
  `TARGET` text,
  `TARGET_BATCH` text,
  `TDATE` datetime DEFAULT NULL,
  `CONC` double DEFAULT NULL,
  `PERCENT_INHIBITION` double DEFAULT NULL,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `PSI` text,
  `CEP_BATCH` text,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `ba_binding_sp_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ba_ca_ic50`
--

DROP TABLE IF EXISTS `ba_ca_ic50`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ba_ca_ic50` (
  `PKEY` text,
  `ABASE_TTORDIRNO` text,
  `ABASE_TTOROBJSEQNO` text,
  `ABASE_TOCCID` text,
  `PROJECT` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `COMPOUND_BATCH` text,
  `TARGET` text,
  `TARGET_BATCH` text,
  `SOURCE` text,
  `TDATE` text,
  `KI` text,
  `DATA_QUALITY` text,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `ba_ca_ic50_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ba_ca_sp`
--

DROP TABLE IF EXISTS `ba_ca_sp`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ba_ca_sp` (
  `PKEY` bigint DEFAULT NULL,
  `ABASE_TTORDIRNO` double DEFAULT NULL,
  `ABASE_TTOROBJSEQNO` double DEFAULT NULL,
  `ABASE_TOCCID` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `PROJECT` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `COMPOUND_BATCH` text,
  `METHOD` text,
  `TARGET` text,
  `TARGET_BATCH` text,
  `TDATE` datetime DEFAULT NULL,
  `CONC` double DEFAULT NULL,
  `PERCENT_CONTROL` double DEFAULT NULL,
  `PERCENT_INHIBITION` double DEFAULT NULL,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `PERCENT_BACKGROUND` double DEFAULT NULL,
  `CREATED_DATE` datetime DEFAULT NULL,
  `ASSAY_TYPE` text,
  `PSI` text,
  `CEP_BATCH` text,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `ba_ca_sp_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ba_camp_efficacy`
--

DROP TABLE IF EXISTS `ba_camp_efficacy`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ba_camp_efficacy` (
  `PKEY` text,
  `ABASE_TTORDIRNO` text,
  `ABASE_TTOROBJSEQNO` text,
  `ABASE_TOCCID` text,
  `PROJECT` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `COMPOUND_BATCH` text,
  `TARGET` text,
  `TARGET_BATCH` text,
  `SOURCE` text,
  `TRACER` text,
  `TDATE` text,
  `KI` text,
  `HILL` text,
  `AGONIST_EFFICACY` text,
  `DATA_QUALITY` text,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `EC50` text,
  `Y_MAX` text,
  `Y_MIN` text,
  `METHOD` text,
  `E_CMAX` text,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `ba_camp_efficacy_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ba_cell_sp`
--

DROP TABLE IF EXISTS `ba_cell_sp`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ba_cell_sp` (
  `PKEY` text,
  `PROJECT` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `TARGET` text,
  `TARGET_BATCH` text,
  `ABASE_TOCCID` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `COMPOUND_BATCH` text,
  `CONC` text,
  `INHIBITION` text,
  `ACT_REL_CONTROL` text,
  `ACT_REL_BACKGROUND` text,
  `COMMENTS` text,
  `OPERATOR` text,
  `TDATE` text,
  `METHOD` text,
  `NOTEBOOK_REF` text,
  `ABASE_TTORDIRNO` text,
  `ABASE_TTOROBJSEQNO` text,
  `ASSAY_TYPE` text,
  `PSI` text,
  `CEP_BATCH` text,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `ba_cell_sp_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ba_enz_inh_comp`
--

DROP TABLE IF EXISTS `ba_enz_inh_comp`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ba_enz_inh_comp` (
  `PKEY` bigint DEFAULT NULL,
  `ABASE_TTORDIRNO` double DEFAULT NULL,
  `ABASE_TTOROBJSEQNO` double DEFAULT NULL,
  `ABASE_TOCCID` text,
  `PROJECT` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `ASSAY_TYPE` text,
  `METHOD` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `COMPOUND_BATCH` text,
  `CMAX` double DEFAULT NULL,
  `TARGET` text,
  `TARGET_BATCH` text,
  `SOURCE` text,
  `SUBSTRATE` text,
  `TDATE` datetime DEFAULT NULL,
  `IC50` double DEFAULT NULL,
  `KM` double DEFAULT NULL,
  `KI` double DEFAULT NULL,
  `HILL` double DEFAULT NULL,
  `Y_MAX` double DEFAULT NULL,
  `Y_MIN` double DEFAULT NULL,
  `E_CMAX` double DEFAULT NULL,
  `DATA_QUALITY` text,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `INHIBITION_TYPE` text,
  `ICFIX50` double DEFAULT NULL,
  `GRAPH` text,
  `ABASE_TTORWELLREFERENCE` text,
  `ABASE_TTORSTATUS` double DEFAULT NULL,
  `EC50` double DEFAULT NULL,
  `PSI` text,
  `CEP_BATCH` text,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `ba_enz_inh_comp_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ba_enzyme_sp`
--

DROP TABLE IF EXISTS `ba_enzyme_sp`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ba_enzyme_sp` (
  `PKEY` bigint DEFAULT NULL,
  `ABASE_TTORDIRNO` double DEFAULT NULL,
  `ABASE_TTOROBJSEQNO` bigint DEFAULT NULL,
  `ABASE_TOCCID` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `PROJECT` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `COMPOUND_BATCH` text,
  `TARGET` text,
  `TARGET_BATCH` text,
  `ASSAY_TYPE` text,
  `TDATE` datetime DEFAULT NULL,
  `SUBSTRATE` text,
  `CONC` double DEFAULT NULL,
  `INHIBITION` double DEFAULT NULL,
  `DATA_QUALITY` text,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `METHOD` double DEFAULT NULL,
  `KM` double DEFAULT NULL,
  `ABASE_TTORWELLREFERENCE` text,
  `INHIBITION_TYPE` text,
  `PSI` text,
  `CEP_BATCH` text,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `ba_enzyme_sp_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ba_functional_dr`
--

DROP TABLE IF EXISTS `ba_functional_dr`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ba_functional_dr` (
  `PKEY` bigint DEFAULT NULL,
  `ABASE_TTORDIRNO` double DEFAULT NULL,
  `ABASE_TTOROBJSEQNO` double DEFAULT NULL,
  `ABASE_TOCCID` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `PROJECT` text,
  `METHOD` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `COMPOUND_BATCH` text,
  `TARGET` text,
  `TARGET_BATCH` text,
  `TDATE` datetime DEFAULT NULL,
  `IC50` double DEFAULT NULL,
  `EC50` double DEFAULT NULL,
  `FKI` double DEFAULT NULL,
  `HILL` double DEFAULT NULL,
  `Y_MAX` double DEFAULT NULL,
  `Y_MIN` double DEFAULT NULL,
  `E_CMAX` double DEFAULT NULL,
  `CMAX` double DEFAULT NULL,
  `DATA_QUALITY` text,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `I_CMAX` double DEFAULT NULL,
  `ASSAY_TYPE` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `GRAPH` text,
  `PSI` text,
  `CEP_BATCH` text,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `ba_functional_dr_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ba_herg_dr`
--

DROP TABLE IF EXISTS `ba_herg_dr`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ba_herg_dr` (
  `PKEY` bigint DEFAULT NULL,
  `ABASE_TTORDIRNO` double DEFAULT NULL,
  `ABASE_TTOROBJSEQNO` bigint DEFAULT NULL,
  `ABASE_TOCCID` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `PROJECT` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `COMPOUND_BATCH` text,
  `TARGET` text,
  `TARGET_BATCH` text,
  `TDATE` datetime DEFAULT NULL,
  `IC50` double DEFAULT NULL,
  `HILL` double DEFAULT NULL,
  `DATA_QUALITY` text,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `GRAPH` text,
  `ASSAY_TYPE` text,
  `R_SQ` double DEFAULT NULL,
  `SB` double DEFAULT NULL,
  `Z` double DEFAULT NULL,
  `REF_CMPD1` double DEFAULT NULL,
  `REF_CMPD2` double DEFAULT NULL,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `ba_herg_dr_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ba_herg_sp`
--

DROP TABLE IF EXISTS `ba_herg_sp`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ba_herg_sp` (
  `PKEY` bigint DEFAULT NULL,
  `ABASE_TTORDIRNO` double DEFAULT NULL,
  `ABASE_TTOROBJSEQNO` bigint DEFAULT NULL,
  `ABASE_TOCCID` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `PROJECT` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `COMPOUND_BATCH` text,
  `TARGET` text,
  `TARGET_BATCH` text,
  `ASSAY_TYPE` text,
  `TDATE` datetime DEFAULT NULL,
  `CONC` double DEFAULT NULL,
  `INHIBITION` double DEFAULT NULL,
  `Z` double DEFAULT NULL,
  `DATA_QUALITY` text,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `ba_herg_sp_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ba_hit_88`
--

DROP TABLE IF EXISTS `ba_hit_88`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ba_hit_88` (
  `PKEY` bigint DEFAULT NULL,
  `ABASE_TTORDIRNO` bigint DEFAULT NULL,
  `ABASE_TTOROBJSEQNO` bigint DEFAULT NULL,
  `PROJECT` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `COMPOUND_BATCH` text,
  `CONC` double DEFAULT NULL,
  `TARGET` text,
  `SOURCE` text,
  `TDATE` datetime DEFAULT NULL,
  `INHIBITION` double DEFAULT NULL,
  `PERCENT_CONTROL` double DEFAULT NULL,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `TARGET_BATCH` text,
  `ABASE_TOCCID` text,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `ba_hit_88_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ba_hts_80_2`
--

DROP TABLE IF EXISTS `ba_hts_80_2`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ba_hts_80_2` (
  `PKEY` text,
  `ABASE_TTORDIRNO` text,
  `ABASE_TTOROBJSEQNO` text,
  `PROJECT` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `COMPOUND_BATCH` text,
  `CONC` text,
  `TARGET` text,
  `SOURCE` text,
  `TDATE` text,
  `RESTACTIVITY` text,
  `SEM` text,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `TARGET_BATCH` text,
  `ABASE_TOCCID` text,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `ba_hts_80_2_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ba_radio_sp`
--

DROP TABLE IF EXISTS `ba_radio_sp`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ba_radio_sp` (
  `PKEY` bigint DEFAULT NULL,
  `ABASE_TTORDIRNO` bigint DEFAULT NULL,
  `ABASE_TTOROBJSEQNO` bigint DEFAULT NULL,
  `ABASE_TTORWELLREFERENCE` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `compound_batch` varchar(16) DEFAULT NULL,
  `CONC` double DEFAULT NULL,
  `INHIBITION` double DEFAULT NULL,
  `HIT` double DEFAULT NULL,
  `COMMENTS` text,
  `ABASE_TOCCID` text,
  `ACTIVATION` double DEFAULT NULL,
  KEY `compound_batch` (`compound_batch`),
  CONSTRAINT `ba_radio_sp_ibfk_1` FOREIGN KEY (`compound_batch`) REFERENCES `bcpvs`.`batch` (`notebook_ref`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ba_rec_ec50`
--

DROP TABLE IF EXISTS `ba_rec_ec50`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ba_rec_ec50` (
  `PKEY` bigint DEFAULT NULL,
  `ABASE_TTORDIRNO` bigint DEFAULT NULL,
  `ABASE_TTOROBJSEQNO` bigint DEFAULT NULL,
  `ABASE_TOCCID` text,
  `PROJECT` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `compound_batch` varchar(16) DEFAULT NULL,
  `METHOD` text,
  `TDATE` datetime DEFAULT NULL,
  `ASSAY_TYPE` text,
  `TARGET` text,
  `TARGET_BATCH` text,
  `SOURCE` text,
  `LIGAND` text,
  `CMAX` double DEFAULT NULL,
  `E_CMAX` double DEFAULT NULL,
  `IC50` double DEFAULT NULL,
  `KD` double DEFAULT NULL,
  `HILL` double DEFAULT NULL,
  `KI` double DEFAULT NULL,
  `Y_MAX` double DEFAULT NULL,
  `Y_MIN` double DEFAULT NULL,
  `DATA_QUALITY` text,
  `OPERATOR` text,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `ICFIX50` double DEFAULT NULL,
  `CREATED_DATE` datetime DEFAULT NULL,
  `GRAPH` text,
  `PSI` text,
  `CEP_BATCH` text,
  KEY `compound_batch` (`compound_batch`),
  CONSTRAINT `ba_rec_ec50_ibfk_1` FOREIGN KEY (`compound_batch`) REFERENCES `bcpvs`.`batch` (`notebook_ref`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `hts_testsets`
--

DROP TABLE IF EXISTS `hts_testsets`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `hts_testsets` (
  `PKEY` bigint DEFAULT NULL,
  `ABASE_TRTSID` text,
  `PROJECT` text,
  `PROTOCOL_ID` text,
  `PROTOCOL_VERSION` text,
  `TDATE` datetime DEFAULT NULL,
  `ASSAY_TYPE` text,
  `TARGET` text,
  `TARGET_BATCH` text,
  `LIGAND` text,
  `KD` text,
  `HIT_THRESHOLD` double DEFAULT NULL,
  `OPERATOR` text,
  `notebook_ref` varchar(16) DEFAULT NULL,
  `COMMENTS` text,
  `METHOD` text
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2022-03-20 13:21:29
