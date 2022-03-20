-- MySQL dump 10.13  Distrib 8.0.26, for Linux (x86_64)
--
-- Host: localhost    Database: bcpvs
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
-- Table structure for table `JCMOL_MOLTABLE`
--

DROP TABLE IF EXISTS `JCMOL_MOLTABLE`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `JCMOL_MOLTABLE` (
  `MOL` longtext CHARACTER SET utf8 COLLATE utf8_bin NOT NULL,
  `COMPOUND_ID` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin NOT NULL,
  UNIQUE KEY `JCMOL_MOLTABLE_ind` (`COMPOUND_ID`),
  CONSTRAINT `JCMOL_MOLTABLE_FK` FOREIGN KEY (`COMPOUND_ID`) REFERENCES `compound` (`compound_id`) ON UPDATE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `JCMOL_MOLTABLE_MOL`
--

DROP TABLE IF EXISTS `JCMOL_MOLTABLE_MOL`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `JCMOL_MOLTABLE_MOL` (
  `COMPOUND_ID` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin NOT NULL,
  `mol` mediumblob,
  UNIQUE KEY `JCMOL_MOLTABLE_MOL_ind` (`COMPOUND_ID`),
  CONSTRAINT `JCMOL_MOLTABLE_MOL_FK` FOREIGN KEY (`COMPOUND_ID`) REFERENCES `JCMOL_MOLTABLE` (`COMPOUND_ID`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `JCMOL_MOLTABLE_MOL_key`
--

DROP TABLE IF EXISTS `JCMOL_MOLTABLE_MOL_key`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `JCMOL_MOLTABLE_MOL_key` (
  `COMPOUND_ID` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin NOT NULL,
  `molkey` binary(192) DEFAULT NULL,
  UNIQUE KEY `JCMOL_MOLTABLE_MOL_key_ind` (`COMPOUND_ID`),
  KEY `JCMOL_MOLTABLE_MOL_key_ind_0` (`molkey`),
  CONSTRAINT `JCMOL_MOLTABLE_MOL_key_FK` FOREIGN KEY (`COMPOUND_ID`) REFERENCES `JCMOL_MOLTABLE_MOL` (`COMPOUND_ID`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `JCMOL_MOLTABLE_MOL_keysim`
--

DROP TABLE IF EXISTS `JCMOL_MOLTABLE_MOL_keysim`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `JCMOL_MOLTABLE_MOL_keysim` (
  `COMPOUND_ID` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin NOT NULL,
  `molkey` binary(192) DEFAULT NULL,
  UNIQUE KEY `JCMOL_MOLTABLE_MOL_keysim_ind` (`COMPOUND_ID`),
  CONSTRAINT `JCMOL_MOLTABLE_MOL_keysim_FK` FOREIGN KEY (`COMPOUND_ID`) REFERENCES `JCMOL_MOLTABLE_MOL` (`COMPOUND_ID`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `batch`
--

DROP TABLE IF EXISTS `batch`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `batch` (
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `notebook_ref` varchar(16) DEFAULT NULL,
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
  `chemspec_regno` varchar(16) DEFAULT NULL,
  `suffix` varchar(20) DEFAULT NULL,
  UNIQUE KEY `batch_idx` (`notebook_ref`),
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `batch_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `batch_invalid`
--

DROP TABLE IF EXISTS `batch_invalid`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `batch_invalid` (
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `notebook_ref` varchar(16) DEFAULT NULL,
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
  `chemspec_regno` varchar(16) DEFAULT NULL,
  `suffix` varchar(20) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `compound`
--

DROP TABLE IF EXISTS `compound`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `compound` (
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `COMPOUND_ID_NUMERIC` bigint DEFAULT NULL,
  `suffix` varchar(5) DEFAULT NULL,
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
  UNIQUE KEY `compound_idx` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `compound_id_sequence`
--

DROP TABLE IF EXISTS `compound_id_sequence`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `compound_id_sequence` (
  `pkey` int NOT NULL,
  PRIMARY KEY (`pkey`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `compound_library`
--

DROP TABLE IF EXISTS `compound_library`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `compound_library` (
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `library_name` varchar(50) DEFAULT NULL,
  `SUPPLIER` text,
  `PROJECT` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `CREATED_BY_PKEY` bigint DEFAULT NULL,
  `UPDATED_DATE` text,
  `UPDATED_BY_PKEY` text,
  `DESCRIPTION` text,
  PRIMARY KEY (`pkey`)
) ENGINE=InnoDB AUTO_INCREMENT=128563351 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `compound_suppliers`
--

DROP TABLE IF EXISTS `compound_suppliers`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `compound_suppliers` (
  `name` varchar(40) DEFAULT NULL,
  `PREFERRED` text
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `compound_type`
--

DROP TABLE IF EXISTS `compound_type`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `compound_type` (
  `type` varchar(20) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `cryst_solvent`
--

DROP TABLE IF EXISTS `cryst_solvent`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `cryst_solvent` (
  `solvent` varchar(20) DEFAULT NULL,
  `MW` double DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `external_id`
--

DROP TABLE IF EXISTS `external_id`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `external_id` (
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `external_id` varchar(100) DEFAULT NULL,
  `CREATED_BY_PKEY` bigint DEFAULT NULL,
  `CREATED_DATE` datetime DEFAULT NULL,
  `UPDATED_BY_PKEY` double DEFAULT NULL,
  `UPDATED_DATE` datetime DEFAULT NULL,
  PRIMARY KEY (`pkey`),
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `external_id_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `external_id_invalid`
--

DROP TABLE IF EXISTS `external_id_invalid`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `external_id_invalid` (
  `PKEY` bigint DEFAULT NULL,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `external_id` varchar(100) DEFAULT NULL,
  `CREATED_BY_PKEY` bigint DEFAULT NULL,
  `CREATED_DATE` datetime DEFAULT NULL,
  `UPDATED_BY_PKEY` double DEFAULT NULL,
  `UPDATED_DATE` datetime DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `id_ok`
--

DROP TABLE IF EXISTS `id_ok`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `id_ok` (
  `ok` varchar(10) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ip_rights_values`
--

DROP TABLE IF EXISTS `ip_rights_values`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ip_rights_values` (
  `ip_rights` varchar(15) DEFAULT NULL,
  `DISP_ORDER` bigint DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `jcsalt_moltable`
--

DROP TABLE IF EXISTS `jcsalt_moltable`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `jcsalt_moltable` (
  `molid` int NOT NULL AUTO_INCREMENT,
  `mol` mediumblob,
  `suffix` varchar(16) CHARACTER SET utf8mb4 COLLATE utf8mb4_0900_ai_ci DEFAULT NULL,
  PRIMARY KEY (`molid`)
) ENGINE=InnoDB AUTO_INCREMENT=33 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `jcsalt_moltable_key`
--

DROP TABLE IF EXISTS `jcsalt_moltable_key`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `jcsalt_moltable_key` (
  `molid` int NOT NULL,
  `molkey` char(192) CHARACTER SET latin1 COLLATE latin1_bin DEFAULT NULL,
  PRIMARY KEY (`molid`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `jcsalt_moltable_keysim`
--

DROP TABLE IF EXISTS `jcsalt_moltable_keysim`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `jcsalt_moltable_keysim` (
  `molid` int NOT NULL,
  `molkey` char(192) CHARACTER SET latin1 COLLATE latin1_bin DEFAULT NULL,
  PRIMARY KEY (`molid`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `jcsalt_moltable_ukey`
--

DROP TABLE IF EXISTS `jcsalt_moltable_ukey`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `jcsalt_moltable_ukey` (
  `molid` int NOT NULL,
  `molkey` mediumblob,
  PRIMARY KEY (`molid`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `lcb_batch`
--

DROP TABLE IF EXISTS `lcb_batch`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `lcb_batch` (
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `notebook_ref` varchar(16) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `library_id_sequence`
--

DROP TABLE IF EXISTS `library_id_sequence`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `library_id_sequence` (
  `pkey` int NOT NULL,
  PRIMARY KEY (`pkey`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `molmass`
--

DROP TABLE IF EXISTS `molmass`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `molmass` (
  `atom` varchar(5) DEFAULT NULL,
  `MASS` double DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `periodic_table`
--

DROP TABLE IF EXISTS `periodic_table`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `periodic_table` (
  `atom` varchar(5) DEFAULT NULL,
  `MASS` double DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `product_type`
--

DROP TABLE IF EXISTS `product_type`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `product_type` (
  `type` varchar(30) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `screening_suppliers`
--

DROP TABLE IF EXISTS `screening_suppliers`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `screening_suppliers` (
  `name` varchar(50) DEFAULT NULL,
  `UPDATEDATE` datetime DEFAULT NULL,
  `COMMENTS` text,
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `PREFERRED` text,
  PRIMARY KEY (`pkey`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `structure_change_log`
--

DROP TABLE IF EXISTS `structure_change_log`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `structure_change_log` (
  `PKEY` bigint DEFAULT NULL,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `NOTEBOOK_REF` text,
  `OLD_COMPOUND_ID` text,
  `NEW_COMPOUND_ID` text,
  `OLD_STRUCTURE` text,
  `NEW_STRUCTURE` text,
  `OLD_PURITY` text,
  `NEW_PURITY` text,
  `OLD_ID_CONFIRMED` text,
  `NEW_ID_CONFIRMED` text,
  `SCHEMA_NAME` text,
  `TABLE_NAME` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `CREATED_BY_PKEY` bigint DEFAULT NULL,
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `structure_change_log_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `compound` (`compound_id`)
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

-- Dump completed on 2022-03-20 13:19:51
