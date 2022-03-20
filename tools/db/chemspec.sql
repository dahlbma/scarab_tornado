-- MySQL dump 10.13  Distrib 8.0.26, for Linux (x86_64)
--
-- Host: localhost    Database: chemspec
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
-- Table structure for table `analytical_sum`
--

DROP TABLE IF EXISTS `analytical_sum`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `analytical_sum` (
  `regno` bigint DEFAULT NULL,
  `FINISHED` text,
  `LOG` text,
  `PKEY` bigint DEFAULT NULL,
  `CREATED_BY_PKEY` bigint DEFAULT NULL,
  `CREATED_DATE` datetime DEFAULT NULL,
  `UPDATED_BY_PKEY` text,
  `UPDATED_DATE` text,
  `RECORD_SOURCE` text,
  `RECORD_SOURCE_PKEY` double DEFAULT NULL,
  KEY `regno` (`regno`),
  CONSTRAINT `analytical_sum_ibfk_1` FOREIGN KEY (`regno`) REFERENCES `chem_info` (`regno`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `chem_info`
--

DROP TABLE IF EXISTS `chem_info`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `chem_info` (
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
  UNIQUE KEY `regno_pk` (`regno`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `chemspec_mol`
--

DROP TABLE IF EXISTS `chemspec_mol`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `chemspec_mol` (
  `molid` int NOT NULL AUTO_INCREMENT,
  `mol` mediumblob,
  `regno` varchar(16) CHARACTER SET utf8mb4 COLLATE utf8mb4_0900_ai_ci DEFAULT NULL,
  PRIMARY KEY (`molid`)
) ENGINE=InnoDB AUTO_INCREMENT=457449 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `chemspec_mol_ukey`
--

DROP TABLE IF EXISTS `chemspec_mol_ukey`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `chemspec_mol_ukey` (
  `molid` int NOT NULL,
  `molkey` mediumblob,
  PRIMARY KEY (`molid`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `chemspec_salt`
--

DROP TABLE IF EXISTS `chemspec_salt`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `chemspec_salt` (
  `suffix` varchar(10) DEFAULT NULL,
  `MOLFILE` text
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `chns_log`
--

DROP TABLE IF EXISTS `chns_log`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `chns_log` (
  `regno` bigint DEFAULT NULL,
  `TDATE` datetime DEFAULT NULL,
  `LOG` text,
  `RAW_DATA_FILE` text,
  `RES_SIGN` text,
  `PROCESSED_DATA_FILE` text,
  `OK` text,
  `EXPERIMENT` text,
  `PURITY` text,
  `PKEY` bigint DEFAULT NULL,
  `METHOD` text,
  `CREATED_BY_PKEY` bigint DEFAULT NULL,
  `CREATED_DATE` datetime DEFAULT NULL,
  `UPDATED_BY_PKEY` text,
  `UPDATED_DATE` text,
  `RECORD_SOURCE` text,
  `RECORD_SOURCE_PKEY` text,
  KEY `regno` (`regno`),
  CONSTRAINT `chns_log_ibfk_1` FOREIGN KEY (`regno`) REFERENCES `chem_info` (`regno`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `color_tbl`
--

DROP TABLE IF EXISTS `color_tbl`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `color_tbl` (
  `color` varchar(20) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `experiment`
--

DROP TABLE IF EXISTS `experiment`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `experiment` (
  `technique` varchar(20) DEFAULT NULL,
  `EXPERIMENT` text,
  `COMMON` double DEFAULT NULL,
  `PURITY_PRIO` double DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ftir_log`
--

DROP TABLE IF EXISTS `ftir_log`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ftir_log` (
  `regno` bigint DEFAULT NULL,
  `TDATE` datetime DEFAULT NULL,
  `LOG` text,
  `RAW_DATA_FILE` text,
  `RES_SIGN` text,
  `PROCESSED_DATA_FILE` text,
  `OK` text,
  `EXPERIMENT` text,
  `PURITY` text,
  `PKEY` bigint DEFAULT NULL,
  `METHOD` text,
  `CREATED_BY_PKEY` bigint DEFAULT NULL,
  `CREATED_DATE` datetime DEFAULT NULL,
  `UPDATED_BY_PKEY` text,
  `UPDATED_DATE` text,
  `RECORD_SOURCE` text,
  `RECORD_SOURCE_PKEY` text,
  KEY `regno` (`regno`),
  CONSTRAINT `ftir_log_ibfk_1` FOREIGN KEY (`regno`) REFERENCES `chem_info` (`regno`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `hplc_uv_log`
--

DROP TABLE IF EXISTS `hplc_uv_log`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `hplc_uv_log` (
  `regno` bigint DEFAULT NULL,
  `TDATE` datetime DEFAULT NULL,
  `LOG` text,
  `RAW_DATA_FILE` text,
  `RES_SIGN` text,
  `PROCESSED_DATA_FILE` text,
  `OK` text,
  `EXPERIMENT` text,
  `PURITY` double DEFAULT NULL,
  `PKEY` bigint DEFAULT NULL,
  `METHOD` text,
  `CREATED_BY_PKEY` bigint DEFAULT NULL,
  `CREATED_DATE` datetime DEFAULT NULL,
  `UPDATED_BY_PKEY` double DEFAULT NULL,
  `UPDATED_DATE` datetime DEFAULT NULL,
  `RECORD_SOURCE` text,
  `RECORD_SOURCE_PKEY` text,
  `RT` text,
  `SYSTEM` text,
  KEY `regno` (`regno`),
  CONSTRAINT `hplc_uv_log_ibfk_1` FOREIGN KEY (`regno`) REFERENCES `chem_info` (`regno`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ms_log`
--

DROP TABLE IF EXISTS `ms_log`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ms_log` (
  `regno` bigint DEFAULT NULL,
  `TDATE` datetime DEFAULT NULL,
  `LOG` text,
  `RAW_DATA_FILE` text,
  `RES_SIGN` text,
  `PROCESSED_DATA_FILE` text,
  `OK` text,
  `EXPERIMENT` text,
  `PURITY` double DEFAULT NULL,
  `PKEY` bigint DEFAULT NULL,
  `EXP_MONO_MASS` double DEFAULT NULL,
  `CALC_MONO_MASS` double DEFAULT NULL,
  `DIFF` double DEFAULT NULL,
  `METHOD` text,
  `METHOD_TEMP` text,
  `CREATED_BY_PKEY` bigint DEFAULT NULL,
  `CREATED_DATE` datetime DEFAULT NULL,
  `UPDATED_BY_PKEY` double DEFAULT NULL,
  `UPDATED_DATE` datetime DEFAULT NULL,
  `INS_SYNC` text,
  `RECORD_SOURCE` text,
  `RECORD_SOURCE_PKEY` text,
  `CONT_REPORT` text,
  KEY `regno` (`regno`),
  CONSTRAINT `ms_log_ibfk_1` FOREIGN KEY (`regno`) REFERENCES `chem_info` (`regno`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `nmr_log`
--

DROP TABLE IF EXISTS `nmr_log`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `nmr_log` (
  `regno` bigint DEFAULT NULL,
  `TDATE` datetime DEFAULT NULL,
  `LOG` text,
  `RAW_DATA_FILE` text,
  `RES_SIGN` text,
  `PROCESSED_DATA_FILE` text,
  `OK` text,
  `EXPERIMENT` text,
  `PURITY` double DEFAULT NULL,
  `PKEY` bigint DEFAULT NULL,
  `METHOD` text,
  `CREATED_BY_PKEY` bigint DEFAULT NULL,
  `CREATED_DATE` datetime DEFAULT NULL,
  `UPDATED_BY_PKEY` double DEFAULT NULL,
  `UPDATED_DATE` datetime DEFAULT NULL,
  `INS_SYNC` text,
  `RECORD_SOURCE` text,
  `RECORD_SOURCE_PKEY` text,
  `CONT_REPORT` text,
  KEY `regno` (`regno`),
  CONSTRAINT `nmr_log_ibfk_1` FOREIGN KEY (`regno`) REFERENCES `chem_info` (`regno`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `other_log`
--

DROP TABLE IF EXISTS `other_log`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `other_log` (
  `regno` bigint DEFAULT NULL,
  `TDATE` datetime DEFAULT NULL,
  `LOG` text,
  `RAW_DATA_FILE` text,
  `PROCESSED_DATA_FILE` text,
  `EXPERIMENT` text,
  `PURITY` text,
  `RES_SIGN` text,
  `OK` text,
  `PKEY` bigint DEFAULT NULL,
  `METHOD` text,
  `CREATED_BY_PKEY` bigint DEFAULT NULL,
  `CREATED_DATE` datetime DEFAULT NULL,
  `UPDATED_BY_PKEY` text,
  `UPDATED_DATE` text,
  `RECORD_SOURCE` text,
  `RECORD_SOURCE_PKEY` text,
  KEY `regno` (`regno`),
  CONSTRAINT `other_log_ibfk_1` FOREIGN KEY (`regno`) REFERENCES `chem_info` (`regno`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `solvent_tbl`
--

DROP TABLE IF EXISTS `solvent_tbl`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `solvent_tbl` (
  `solvent` varchar(20) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `spec_analysis`
--

DROP TABLE IF EXISTS `spec_analysis`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `spec_analysis` (
  `regno` bigint DEFAULT NULL,
  `RECEIVED_BY` text,
  `RECEIVED_DATE` datetime DEFAULT NULL,
  `ID` text,
  `APPROVED_BY` text,
  `APPROVED_DATE` text,
  `BINDER` text,
  `PKEY` bigint DEFAULT NULL,
  `PLATE_ID` text,
  `CREATED_BY_PKEY` bigint DEFAULT NULL,
  `CREATED_DATE` datetime DEFAULT NULL,
  `UPDATED_BY_PKEY` text,
  `UPDATED_DATE` text,
  `SELF_PKEY` text,
  `RECORD_SOURCE` text,
  `RECORD_SOURCE_PKEY` bigint DEFAULT NULL,
  KEY `regno` (`regno`),
  CONSTRAINT `spec_analysis_ibfk_1` FOREIGN KEY (`regno`) REFERENCES `chem_info` (`regno`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `storage_tbl`
--

DROP TABLE IF EXISTS `storage_tbl`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `storage_tbl` (
  `storage` varchar(30) DEFAULT NULL
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

-- Dump completed on 2022-03-20 13:21:15
