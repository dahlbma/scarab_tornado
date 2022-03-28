-- MySQL dump 10.13  Distrib 8.0.26, for Linux (x86_64)
--
-- Host: localhost    Database: assay
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
-- Table structure for table `assay_types`
--

DROP TABLE IF EXISTS `assay_types`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `assay_types` (
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `assay_type` varchar(100) DEFAULT NULL,
  PRIMARY KEY (`pkey`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ba_bf_camp_dr`
--

DROP TABLE IF EXISTS `ba_bf_camp_dr`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ba_bf_camp_dr` (
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `project` varchar(40) DEFAULT NULL,
  `method` varchar(20) DEFAULT NULL,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `compound_batch` varchar(16) DEFAULT NULL,
  `target` varchar(20) DEFAULT NULL,
  `TARGET_BATCH` text,
  `TDATE` datetime DEFAULT NULL,
  `IC50` text,
  `EC50` double DEFAULT NULL,
  `HILL` double DEFAULT NULL,
  `Y_MAX` double DEFAULT NULL,
  `M_MIN` double DEFAULT NULL,
  `DATA_QUALITY` text,
  `operator` varchar(20) DEFAULT NULL,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `FKI` text,
  PRIMARY KEY (`pkey`),
  KEY `compound_batch` (`compound_batch`),
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `ba_bf_camp_dr_ibfk_1` FOREIGN KEY (`compound_batch`) REFERENCES `bcpvs`.`batch` (`notebook_ref`),
  CONSTRAINT `ba_bf_camp_dr_ibfk_2` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ba_bf_camp_sp`
--

DROP TABLE IF EXISTS `ba_bf_camp_sp`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ba_bf_camp_sp` (
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `project` varchar(20) DEFAULT NULL,
  `method` varchar(20) DEFAULT NULL,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `compound_batch` varchar(16) DEFAULT NULL,
  `target` varchar(20) DEFAULT NULL,
  `TARGET_BATCH` text,
  `TDATE` datetime DEFAULT NULL,
  `CONC` double DEFAULT NULL,
  `PERCENT_INHIBITION` text,
  `PERCENT_BACKGROUND` double DEFAULT NULL,
  `operator` varchar(20) DEFAULT NULL,
  `notebook_ref` varchar(20) DEFAULT NULL,
  `COMMENTS` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `PERCENT_CONTROL` text,
  PRIMARY KEY (`pkey`),
  KEY `compound_batch` (`compound_batch`),
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `ba_bf_camp_sp_ibfk_1` FOREIGN KEY (`compound_batch`) REFERENCES `bcpvs`.`batch` (`notebook_ref`),
  CONSTRAINT `ba_bf_camp_sp_ibfk_2` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ba_biofocus_dr`
--

DROP TABLE IF EXISTS `ba_biofocus_dr`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ba_biofocus_dr` (
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `project` varchar(20) DEFAULT NULL,
  `method` varchar(20) DEFAULT NULL,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `compound_batch` varchar(16) DEFAULT NULL,
  `target` varchar(20) DEFAULT NULL,
  `TDATE` datetime DEFAULT NULL,
  `IC50` double DEFAULT NULL,
  `HILL` double DEFAULT NULL,
  `Y_MAX` double DEFAULT NULL,
  `M_MIN` double DEFAULT NULL,
  `DATA_QUALITY` text,
  `operator` varchar(20) DEFAULT NULL,
  `notebook_ref` varchar(20) DEFAULT NULL,
  `COMMENTS` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `EC50` double DEFAULT NULL,
  `TARGET_BATCH` text,
  `FKI` text,
  PRIMARY KEY (`pkey`),
  KEY `compound_batch` (`compound_batch`),
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `ba_biofocus_dr_ibfk_1` FOREIGN KEY (`compound_batch`) REFERENCES `bcpvs`.`batch` (`notebook_ref`),
  CONSTRAINT `ba_biofocus_dr_ibfk_2` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ba_biofocus_sp`
--

DROP TABLE IF EXISTS `ba_biofocus_sp`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ba_biofocus_sp` (
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `project` varchar(20) DEFAULT NULL,
  `method` varchar(20) DEFAULT NULL,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `compound_batch` varchar(16) DEFAULT NULL,
  `target` varchar(20) DEFAULT NULL,
  `TDATE` datetime DEFAULT NULL,
  `CONC` double DEFAULT NULL,
  `PERCENT_INHIBITION` double DEFAULT NULL,
  `PERCENT_BACKGROUND` double DEFAULT NULL,
  `operator` varchar(20) DEFAULT NULL,
  `notebook_ref` varchar(20) DEFAULT NULL,
  `COMMENTS` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `TARGET_BATCH` text,
  `PERCENT_CONTROL` double DEFAULT NULL,
  PRIMARY KEY (`pkey`),
  KEY `compound_batch` (`compound_batch`),
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `ba_biofocus_sp_ibfk_1` FOREIGN KEY (`compound_batch`) REFERENCES `bcpvs`.`batch` (`notebook_ref`),
  CONSTRAINT `ba_biofocus_sp_ibfk_2` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ba_mds_dr`
--

DROP TABLE IF EXISTS `ba_mds_dr`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ba_mds_dr` (
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `compound_batch` varchar(16) DEFAULT NULL,
  `project` varchar(20) DEFAULT NULL,
  `CRO_COMPOUND_ID` text,
  `EXPERIMENT_DATE` datetime DEFAULT NULL,
  `CREATED_DATE` datetime DEFAULT NULL,
  `DATE_OF_STUDY_END` datetime DEFAULT NULL,
  `CRO_SUDY_ID` bigint DEFAULT NULL,
  `MOLWEIGHT` double DEFAULT NULL,
  `CAT_NO` text,
  `target` varchar(80) DEFAULT NULL,
  `LIGAND` text,
  `species` varchar(20) DEFAULT NULL,
  `SOURCE` text,
  `VEHICLE` text,
  `CONCENTRATION` bigint DEFAULT NULL,
  `INHIBITION` bigint DEFAULT NULL,
  `IC50` text,
  `KI` text,
  `NH` text,
  `AGONIST_ACTIVATION_PCT` text,
  `ANTAGONIST_INHIBITION_PCT` double DEFAULT NULL,
  `METHOD` text,
  `operator` varchar(20) DEFAULT NULL,
  `notebook_ref` varchar(20) DEFAULT NULL,
  `CRO_REPORT` text,
  `BWISE_REPORT` text,
  `COMMENTS` text,
  PRIMARY KEY (`pkey`),
  KEY `compound_batch` (`compound_batch`),
  CONSTRAINT `ba_mds_dr_ibfk_1` FOREIGN KEY (`compound_batch`) REFERENCES `bcpvs`.`batch` (`notebook_ref`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `cerep_assays`
--

DROP TABLE IF EXISTS `cerep_assays`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `cerep_assays` (
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `ASSAY` text,
  `REF_COMPOUND` text,
  `REF_IC50` text,
  `REF_EC50` text,
  `COMMENTS` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `CREATED_BY_PKEY` bigint DEFAULT NULL,
  `UPDATED_DATE` text,
  `UPDATED_BY_PKEY` text,
  PRIMARY KEY (`pkey`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `cerep_functional_dr`
--

DROP TABLE IF EXISTS `cerep_functional_dr`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `cerep_functional_dr` (
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `project` varchar(20) DEFAULT NULL,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `compound_batch` varchar(16) DEFAULT NULL,
  `target` varchar(80) DEFAULT NULL,
  `TARGET_BATCH` text,
  `TDATE` datetime DEFAULT NULL,
  `IC50` double DEFAULT NULL,
  `EC50` double DEFAULT NULL,
  `GRAPH` text,
  `operator` varchar(20) DEFAULT NULL,
  `notebook_ref` varchar(20) DEFAULT NULL,
  `COMMENTS` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `CREATED_BY_PKEY` bigint DEFAULT NULL,
  `UPDATED_DATE` text,
  `UPDATED_BY_PKEY` text,
  PRIMARY KEY (`pkey`),
  KEY `compound_batch` (`compound_batch`),
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `cerep_functional_dr_ibfk_1` FOREIGN KEY (`compound_batch`) REFERENCES `bcpvs`.`batch` (`notebook_ref`),
  CONSTRAINT `cerep_functional_dr_ibfk_2` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `cerep_functional_sp`
--

DROP TABLE IF EXISTS `cerep_functional_sp`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `cerep_functional_sp` (
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `project` varchar(20) DEFAULT NULL,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `compound_batch` varchar(16) DEFAULT NULL,
  `target` varchar(80) DEFAULT NULL,
  `TDATE` datetime DEFAULT NULL,
  `CONC` double DEFAULT NULL,
  `N` double DEFAULT NULL,
  `INHIBITION` double DEFAULT NULL,
  `STIMULATION` double DEFAULT NULL,
  `operator` varchar(20) DEFAULT NULL,
  `notebook_ref` varchar(20) DEFAULT NULL,
  `COMMENTS` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `CREATED_BY_PKEY` bigint DEFAULT NULL,
  `UPDATED_DATE` text,
  `UPDATED_BY_PKEY` text,
  `TARGET_BATCH` text,
  PRIMARY KEY (`pkey`),
  KEY `compound_batch` (`compound_batch`),
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `cerep_functional_sp_ibfk_1` FOREIGN KEY (`compound_batch`) REFERENCES `bcpvs`.`batch` (`notebook_ref`),
  CONSTRAINT `cerep_functional_sp_ibfk_2` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `cerep_ki`
--

DROP TABLE IF EXISTS `cerep_ki`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `cerep_ki` (
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `project` varchar(20) DEFAULT NULL,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `compound_batch` varchar(16) DEFAULT NULL,
  `target` varchar(80) DEFAULT NULL,
  `TDATE` datetime DEFAULT NULL,
  `KI` double DEFAULT NULL,
  `operator` varchar(20) DEFAULT NULL,
  `notebook_ref` varchar(20) DEFAULT NULL,
  `COMMENTS` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `CREATED_BY_PKEY` bigint DEFAULT NULL,
  `UPDATED_DATE` text,
  `UPDATED_BY_PKEY` text,
  `TARGET_BATCH` text,
  `GRAPH` text,
  PRIMARY KEY (`pkey`),
  KEY `compound_batch` (`compound_batch`),
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `cerep_ki_ibfk_1` FOREIGN KEY (`compound_batch`) REFERENCES `bcpvs`.`batch` (`notebook_ref`),
  CONSTRAINT `cerep_ki_ibfk_2` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `cerep_receptors`
--

DROP TABLE IF EXISTS `cerep_receptors`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `cerep_receptors` (
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `RECEPTOR` text,
  `ORIGIN` text,
  `REF_COMPOUND` text,
  `LIT_REF` text,
  `LIGAND` text,
  `LIGAND_CONC` text,
  `NONSPEC_COMPOUND` text,
  `INCUBATION` text,
  `COMMENTS` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `CREATED_BY_PKEY` bigint DEFAULT NULL,
  `UPDATED_DATE` text,
  `UPDATED_BY_PKEY` text,
  PRIMARY KEY (`pkey`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `cerep_screen`
--

DROP TABLE IF EXISTS `cerep_screen`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `cerep_screen` (
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `project` varchar(20) DEFAULT NULL,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `compound_batch` varchar(16) DEFAULT NULL,
  `target` varchar(80) DEFAULT NULL,
  `TDATE` datetime DEFAULT NULL,
  `CONC` double DEFAULT NULL,
  `N` bigint DEFAULT NULL,
  `INHIBITION` double DEFAULT NULL,
  `operator` varchar(20) DEFAULT NULL,
  `notebook_ref` varchar(20) DEFAULT NULL,
  `COMMENTS` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `CREATED_BY_PKEY` bigint DEFAULT NULL,
  `UPDATED_DATE` text,
  `UPDATED_BY_PKEY` text,
  `TARGET_BATCH` text,
  PRIMARY KEY (`pkey`),
  KEY `compound_batch` (`compound_batch`),
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `cerep_screen_ibfk_1` FOREIGN KEY (`compound_batch`) REFERENCES `bcpvs`.`batch` (`notebook_ref`),
  CONSTRAINT `cerep_screen_ibfk_2` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `detection_types`
--

DROP TABLE IF EXISTS `detection_types`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `detection_types` (
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `DETECTION_TYPE` text,
  PRIMARY KEY (`pkey`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `lcb_dr`
--

DROP TABLE IF EXISTS `lcb_dr`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `lcb_dr` (
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `project` varchar(40) DEFAULT NULL,
  `TARGET` text,
  `PLATE_ID` text,
  `compound_batch` varchar(16) DEFAULT NULL,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `assay_type` varchar(20) DEFAULT NULL,
  `detection_type` varchar(20) DEFAULT NULL,
  `CMAX` double DEFAULT NULL,
  `Y_MAX` double DEFAULT NULL,
  `M_MIN` double DEFAULT NULL,
  `HILL` double DEFAULT NULL,
  `IC50` double DEFAULT NULL,
  `EC50` double DEFAULT NULL,
  `I_CMAX` double DEFAULT NULL,
  `E_CMAX` double DEFAULT NULL,
  `GRAPH` text,
  `TDATE` datetime DEFAULT NULL,
  `operator` varchar(20) DEFAULT NULL,
  `eln_id` varchar(20) DEFAULT NULL,
  `COMMENTS` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `CONFIRMED` text,
  PRIMARY KEY (`pkey`),
  KEY `compound_batch` (`compound_batch`),
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `lcb_dr_ibfk_1` FOREIGN KEY (`compound_batch`) REFERENCES `bcpvs`.`batch` (`notebook_ref`),
  CONSTRAINT `lcb_dr_ibfk_2` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `lcb_dr_invalid`
--

DROP TABLE IF EXISTS `lcb_dr_invalid`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `lcb_dr_invalid` (
  `PKEY` bigint DEFAULT NULL,
  `project` varchar(40) DEFAULT NULL,
  `TARGET` text,
  `PLATE_ID` text,
  `compound_batch` varchar(16) DEFAULT NULL,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `assay_type` varchar(20) DEFAULT NULL,
  `detection_type` varchar(20) DEFAULT NULL,
  `CMAX` double DEFAULT NULL,
  `Y_MAX` double DEFAULT NULL,
  `M_MIN` double DEFAULT NULL,
  `HILL` double DEFAULT NULL,
  `IC50` double DEFAULT NULL,
  `EC50` double DEFAULT NULL,
  `I_CMAX` double DEFAULT NULL,
  `E_CMAX` double DEFAULT NULL,
  `GRAPH` text,
  `TDATE` datetime DEFAULT NULL,
  `operator` varchar(20) DEFAULT NULL,
  `eln_id` varchar(20) DEFAULT NULL,
  `COMMENTS` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `CONFIRMED` text
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `lcb_sp`
--

DROP TABLE IF EXISTS `lcb_sp`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `lcb_sp` (
  `PKEY` bigint DEFAULT NULL,
  `project` varchar(40) DEFAULT NULL,
  `target` varchar(40) DEFAULT NULL,
  `PLATE_ID` text,
  `compound_batch` varchar(16) DEFAULT NULL,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `assay_type` varchar(20) DEFAULT NULL,
  `detection_type` varchar(20) DEFAULT NULL,
  `CONC` double DEFAULT NULL,
  `INHIBITION` double DEFAULT NULL,
  `ACTIVATION` double DEFAULT NULL,
  `HIT` double DEFAULT NULL,
  `HIT_THRESHOLD` double DEFAULT NULL,
  `TDATE` datetime DEFAULT NULL,
  `operator` varchar(40) DEFAULT NULL,
  `eln_id` varchar(40) DEFAULT NULL,
  `COMMENTS` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `WELL_ID` text,
  KEY `compound_batch` (`compound_batch`),
  CONSTRAINT `lcb_sp_ibfk_1` FOREIGN KEY (`compound_batch`) REFERENCES `bcpvs`.`batch` (`notebook_ref`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `lcb_sp_invalid`
--

DROP TABLE IF EXISTS `lcb_sp_invalid`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `lcb_sp_invalid` (
  `PKEY` bigint DEFAULT NULL,
  `project` varchar(40) DEFAULT NULL,
  `target` varchar(40) DEFAULT NULL,
  `PLATE_ID` text,
  `compound_batch` varchar(16) DEFAULT NULL,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `assay_type` varchar(20) DEFAULT NULL,
  `detection_type` varchar(20) DEFAULT NULL,
  `CONC` double DEFAULT NULL,
  `INHIBITION` double DEFAULT NULL,
  `ACTIVATION` double DEFAULT NULL,
  `HIT` double DEFAULT NULL,
  `HIT_THRESHOLD` double DEFAULT NULL,
  `TDATE` datetime DEFAULT NULL,
  `operator` varchar(40) DEFAULT NULL,
  `eln_id` varchar(40) DEFAULT NULL,
  `COMMENTS` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `WELL_ID` text
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `lcb_sp_schneider_cysm`
--

DROP TABLE IF EXISTS `lcb_sp_schneider_cysm`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `lcb_sp_schneider_cysm` (
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `project` varchar(40) DEFAULT NULL,
  `target` varchar(80) DEFAULT NULL,
  `PLATE_ID` text,
  `compound_batch` varchar(16) DEFAULT NULL,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `assay_type` varchar(20) DEFAULT NULL,
  `detection_type` varchar(20) DEFAULT NULL,
  `CONC` double DEFAULT NULL,
  `INHIBITION` text,
  `ACTIVATION` double DEFAULT NULL,
  `HIT` bigint DEFAULT NULL,
  `HIT_THRESHOLD` double DEFAULT NULL,
  `TDATE` datetime DEFAULT NULL,
  `operator` varchar(20) DEFAULT NULL,
  `eln_id` varchar(20) DEFAULT NULL,
  `COMMENTS` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `WELL_ID` text,
  PRIMARY KEY (`pkey`),
  KEY `compound_batch` (`compound_batch`),
  CONSTRAINT `lcb_sp_schneider_cysm_ibfk_1` FOREIGN KEY (`compound_batch`) REFERENCES `bcpvs`.`batch` (`notebook_ref`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `lcb_sp_schneider_cysm_invalid`
--

DROP TABLE IF EXISTS `lcb_sp_schneider_cysm_invalid`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `lcb_sp_schneider_cysm_invalid` (
  `PKEY` bigint DEFAULT NULL,
  `project` varchar(40) DEFAULT NULL,
  `target` varchar(80) DEFAULT NULL,
  `PLATE_ID` text,
  `compound_batch` varchar(16) DEFAULT NULL,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `assay_type` varchar(20) DEFAULT NULL,
  `detection_type` varchar(20) DEFAULT NULL,
  `CONC` double DEFAULT NULL,
  `INHIBITION` text,
  `ACTIVATION` double DEFAULT NULL,
  `HIT` bigint DEFAULT NULL,
  `HIT_THRESHOLD` double DEFAULT NULL,
  `TDATE` datetime DEFAULT NULL,
  `operator` varchar(20) DEFAULT NULL,
  `eln_id` varchar(20) DEFAULT NULL,
  `COMMENTS` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `WELL_ID` text
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `metabolite`
--

DROP TABLE IF EXISTS `metabolite`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `metabolite` (
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `project` varchar(40) DEFAULT NULL,
  `TDATE` datetime DEFAULT NULL,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `compound_batch` varchar(16) DEFAULT NULL,
  `EXP_METHOD` double DEFAULT NULL,
  `INC_METHOD` double DEFAULT NULL,
  `RESULT_DOC` bigint DEFAULT NULL,
  `DATA_QUALITY` text,
  `operator` varchar(20) DEFAULT NULL,
  `notebook_ref` varchar(20) DEFAULT NULL,
  `CREATED_DATE` datetime DEFAULT NULL,
  `CREATED_BY_PKEY` bigint DEFAULT NULL,
  `UPDATED_DATE` text,
  `UPDATED_BY_PKEY` text,
  `EXPERIMENT_TYPE` text,
  `BINDER_NUMBER` text,
  `COMMENTS` text,
  PRIMARY KEY (`pkey`),
  KEY `compound_batch` (`compound_batch`),
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `metabolite_ibfk_1` FOREIGN KEY (`compound_batch`) REFERENCES `bcpvs`.`batch` (`notebook_ref`),
  CONSTRAINT `metabolite_ibfk_2` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `upstate_sp`
--

DROP TABLE IF EXISTS `upstate_sp`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `upstate_sp` (
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `project` varchar(40) DEFAULT NULL,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `compound_batch` varchar(16) DEFAULT NULL,
  `target` varchar(80) DEFAULT NULL,
  `method` varchar(20) DEFAULT NULL,
  `TDATE` datetime DEFAULT NULL,
  `CONC` double DEFAULT NULL,
  `PERCENT_CONTROL` double DEFAULT NULL,
  `DATA_QUALITY` text,
  `operator` varchar(20) DEFAULT NULL,
  `NOTEBOOK_REF` text,
  `COMMENTS` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `CREATED_BY_PKEY` bigint DEFAULT NULL,
  `UPDATED_DATE` text,
  `UPDATED_BY_PKEY` text,
  `SUBSTRATE_CONC` double DEFAULT NULL,
  PRIMARY KEY (`pkey`),
  KEY `compound_batch` (`compound_batch`),
  KEY `compound_id` (`compound_id`),
  CONSTRAINT `upstate_sp_ibfk_1` FOREIGN KEY (`compound_batch`) REFERENCES `bcpvs`.`batch` (`notebook_ref`),
  CONSTRAINT `upstate_sp_ibfk_2` FOREIGN KEY (`compound_id`) REFERENCES `bcpvs`.`compound` (`compound_id`)
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

-- Dump completed on 2022-03-20 13:20:04
