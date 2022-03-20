-- MySQL dump 10.13  Distrib 8.0.26, for Linux (x86_64)
--
-- Host: localhost    Database: chem_reg
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
-- Table structure for table `CHEM`
--

DROP TABLE IF EXISTS `CHEM`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `CHEM` (
  `MOL` longtext CHARACTER SET utf8 COLLATE utf8_bin NOT NULL,
  `REGNO` int NOT NULL,
  UNIQUE KEY `CHEM_ind` (`REGNO`),
  CONSTRAINT `CHEM_ibfk_1` FOREIGN KEY (`REGNO`) REFERENCES `chem_info` (`regno`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `CHEM_MOL`
--

DROP TABLE IF EXISTS `CHEM_MOL`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `CHEM_MOL` (
  `REGNO` int NOT NULL,
  `mol` mediumblob,
  UNIQUE KEY `CHEM_MOL_ind` (`REGNO`),
  CONSTRAINT `CHEM_MOL_FK` FOREIGN KEY (`REGNO`) REFERENCES `CHEM` (`REGNO`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `CHEM_MOL_key`
--

DROP TABLE IF EXISTS `CHEM_MOL_key`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `CHEM_MOL_key` (
  `REGNO` int NOT NULL,
  `molkey` binary(192) DEFAULT NULL,
  UNIQUE KEY `CHEM_MOL_key_ind` (`REGNO`),
  KEY `CHEM_MOL_key_ind_0` (`molkey`),
  CONSTRAINT `CHEM_MOL_key_FK` FOREIGN KEY (`REGNO`) REFERENCES `CHEM_MOL` (`REGNO`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `CHEM_MOL_keysim`
--

DROP TABLE IF EXISTS `CHEM_MOL_keysim`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `CHEM_MOL_keysim` (
  `REGNO` int NOT NULL,
  `molkey` binary(192) DEFAULT NULL,
  UNIQUE KEY `CHEM_MOL_keysim_ind` (`REGNO`),
  CONSTRAINT `CHEM_MOL_keysim_FK` FOREIGN KEY (`REGNO`) REFERENCES `CHEM_MOL` (`REGNO`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `chem_info`
--

DROP TABLE IF EXISTS `chem_info`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `chem_info` (
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
  `suffix` varchar(20) DEFAULT NULL,
  `IP_RIGHTS` varchar(20) DEFAULT NULL,
  `sdfile_sequence` int DEFAULT NULL,
  UNIQUE KEY `regno_pk` (`regno`),
  UNIQUE KEY `index_name` (`JPAGE`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `chem_info_mol`
--

DROP TABLE IF EXISTS `chem_info_mol`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `chem_info_mol` (
  `molid` int NOT NULL AUTO_INCREMENT,
  `mol` mediumblob,
  `regno` varchar(16) CHARACTER SET utf8mb4 COLLATE utf8mb4_0900_ai_ci DEFAULT NULL,
  PRIMARY KEY (`molid`),
  KEY `regno` (`regno`)
) ENGINE=InnoDB AUTO_INCREMENT=42249 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `chem_info_mol_key`
--

DROP TABLE IF EXISTS `chem_info_mol_key`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `chem_info_mol_key` (
  `molid` int NOT NULL,
  `molkey` char(192) DEFAULT NULL,
  PRIMARY KEY (`molid`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `chem_info_mol_ukey`
--

DROP TABLE IF EXISTS `chem_info_mol_ukey`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `chem_info_mol_ukey` (
  `molid` int NOT NULL,
  `molkey` mediumblob,
  PRIMARY KEY (`molid`),
  KEY `mol_ukey_uidx` (`molkey`(255))
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `chemreg_dist`
--

DROP TABLE IF EXISTS `chemreg_dist`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `chemreg_dist` (
  `pkey` int NOT NULL AUTO_INCREMENT,
  `program` longblob,
  `os` varchar(64) DEFAULT NULL,
  `version` varchar(64) DEFAULT NULL,
  PRIMARY KEY (`pkey`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `regno_sequence`
--

DROP TABLE IF EXISTS `regno_sequence`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `regno_sequence` (
  `pkey` int NOT NULL,
  PRIMARY KEY (`pkey`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `salts`
--

DROP TABLE IF EXISTS `salts`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `salts` (
  `pkey` int NOT NULL AUTO_INCREMENT,
  `suffix` varchar(10) DEFAULT NULL,
  `MOLFILE` mediumtext,
  `smiles` varchar(200) DEFAULT NULL,
  `mf` varchar(200) DEFAULT NULL,
  `mw` float DEFAULT NULL,
  PRIMARY KEY (`pkey`)
) ENGINE=InnoDB AUTO_INCREMENT=82 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `salts_sequence`
--

DROP TABLE IF EXISTS `salts_sequence`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `salts_sequence` (
  `pkey` int NOT NULL,
  PRIMARY KEY (`pkey`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `sdfile_sequence`
--

DROP TABLE IF EXISTS `sdfile_sequence`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `sdfile_sequence` (
  `pkey` int NOT NULL,
  PRIMARY KEY (`pkey`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tmp_mol`
--

DROP TABLE IF EXISTS `tmp_mol`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `tmp_mol` (
  `pkey` int DEFAULT NULL,
  `molfile` text
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tmp_mol_sequence`
--

DROP TABLE IF EXISTS `tmp_mol_sequence`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `tmp_mol_sequence` (
  `pkey` int NOT NULL,
  PRIMARY KEY (`pkey`)
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

-- Dump completed on 2022-03-20 20:39:09
