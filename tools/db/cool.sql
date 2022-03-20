-- MySQL dump 10.13  Distrib 8.0.26, for Linux (x86_64)
--
-- Host: localhost    Database: cool
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
-- Temporary view structure for view `V_PICK_RECIPIENT`
--

DROP TABLE IF EXISTS `V_PICK_RECIPIENT`;
/*!50001 DROP VIEW IF EXISTS `V_PICK_RECIPIENT`*/;
SET @saved_cs_client     = @@character_set_client;
/*!50503 SET character_set_client = utf8mb4 */;
/*!50001 CREATE VIEW `V_PICK_RECIPIENT` AS SELECT 
 1 AS `PKEY`,
 1 AS `NAME`*/;
SET character_set_client = @saved_cs_client;

--
-- Temporary view structure for view `V_PLATE_RNAI`
--

DROP TABLE IF EXISTS `V_PLATE_RNAI`;
/*!50001 DROP VIEW IF EXISTS `V_PLATE_RNAI`*/;
SET @saved_cs_client     = @@character_set_client;
/*!50503 SET character_set_client = utf8mb4 */;
/*!50001 CREATE VIEW `V_PLATE_RNAI` AS SELECT 
 1 AS `PLATE_ID`,
 1 AS `WELL`,
 1 AS `VENDOR_REAGENT_ID`,
 1 AS `CONC`*/;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `config`
--

DROP TABLE IF EXISTS `config`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `config` (
  `config_id` varchar(7) DEFAULT NULL,
  `well` varchar(10) DEFAULT NULL,
  `compound_id` varchar(26) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `notebook_ref` varchar(42) DEFAULT NULL,
  `FORM` text,
  `CONC` double DEFAULT NULL,
  `VOLUME` int DEFAULT NULL,
  KEY `idx_configid` (`config_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `config_key`
--

DROP TABLE IF EXISTS `config_key`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `config_key` (
  `config_id` varchar(7) NOT NULL,
  `KEY` text,
  PRIMARY KEY (`config_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `map96to384`
--

DROP TABLE IF EXISTS `map96to384`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `map96to384` (
  `QUADRANT` bigint DEFAULT NULL,
  `well96` varchar(3) DEFAULT NULL,
  `well384` varchar(4) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `plate`
--

DROP TABLE IF EXISTS `plate`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `plate` (
  `plate_id` varchar(7) NOT NULL,
  `config_id` varchar(7) DEFAULT NULL,
  `TYPE_ID` bigint DEFAULT NULL,
  `DILUTION` bigint DEFAULT NULL,
  `location` varchar(100) DEFAULT NULL,
  `COMMENTS` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `CREATOR_PKEY` bigint DEFAULT NULL,
  `RECIPIENT_PKEY` double DEFAULT NULL,
  `UPDATED_DATE` datetime DEFAULT NULL,
  `UPDATER_PKEY` bigint DEFAULT NULL,
  `DISCARDED` varchar(10) DEFAULT NULL,
  PRIMARY KEY (`plate_id`),
  KEY `idx_plate_config_id` (`config_id`),
  KEY `idx_plate_plate_id` (`plate_id`),
  KEY `idx_plate_TYPE_ID` (`TYPE_ID`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `plate_sequence`
--

DROP TABLE IF EXISTS `plate_sequence`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `plate_sequence` (
  `pkey` int NOT NULL,
  PRIMARY KEY (`pkey`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `plate_type`
--

DROP TABLE IF EXISTS `plate_type`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `plate_type` (
  `type_id` varchar(7) NOT NULL,
  `NAME` text,
  `WELLS` bigint DEFAULT NULL,
  `TERMINATED_DATE` text,
  `PLFTID` text,
  `ALLOW_DRAW_PLATE_ID` double DEFAULT NULL,
  `LABEL_FORMAT` text,
  `RACKS_PER_PLATE` double DEFAULT NULL,
  PRIMARY KEY (`type_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `plating_sequence`
--

DROP TABLE IF EXISTS `plating_sequence`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `plating_sequence` (
  `TYPE_ID` bigint DEFAULT NULL,
  `SEQ` bigint DEFAULT NULL,
  `well` varchar(10) DEFAULT NULL,
  KEY `idx_plating_sequence_SEQ` (`SEQ`),
  KEY `idx_plating_sequence_TYPE_ID` (`TYPE_ID`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `solubility_problem`
--

DROP TABLE IF EXISTS `solubility_problem`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `solubility_problem` (
  `notebook_ref` varchar(20) DEFAULT NULL,
  `TESTED_CONC` double DEFAULT NULL,
  `TDATE` datetime DEFAULT NULL,
  `OPERATOR` text,
  `pkey` bigint NOT NULL,
  `PROBLEM` text,
  PRIMARY KEY (`pkey`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Temporary view structure for view `v_plate`
--

DROP TABLE IF EXISTS `v_plate`;
/*!50001 DROP VIEW IF EXISTS `v_plate`*/;
SET @saved_cs_client     = @@character_set_client;
/*!50503 SET character_set_client = utf8mb4 */;
/*!50001 CREATE VIEW `v_plate` AS SELECT 
 1 AS `PLATE_ID`,
 1 AS `WELL`,
 1 AS `COMPOUND_ID`,
 1 AS `NOTEBOOK_REF`,
 1 AS `FORM`,
 1 AS `CONC`*/;
SET character_set_client = @saved_cs_client;

--
-- Final view structure for view `V_PICK_RECIPIENT`
--

/*!50001 DROP VIEW IF EXISTS `V_PICK_RECIPIENT`*/;
/*!50001 SET @saved_cs_client          = @@character_set_client */;
/*!50001 SET @saved_cs_results         = @@character_set_results */;
/*!50001 SET @saved_col_connection     = @@collation_connection */;
/*!50001 SET character_set_client      = utf8mb4 */;
/*!50001 SET character_set_results     = utf8mb4 */;
/*!50001 SET collation_connection      = utf8mb4_0900_ai_ci */;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`root`@`localhost` SQL SECURITY DEFINER */
/*!50001 VIEW `V_PICK_RECIPIENT` (`PKEY`,`NAME`) AS select `hive`.`user_details`.`pkey` AS `pkey`,`hive`.`user_details`.`fullname` AS `fullname` from `hive`.`user_details` where ((`hive`.`user_details`.`HIVELOCATION` is not null) and ((`hive`.`user_details`.`TERMINATED_DATE` is null) or (`hive`.`user_details`.`TERMINATED_DATE` > curdate()))) order by lower(`hive`.`user_details`.`fullname`) */;
/*!50001 SET character_set_client      = @saved_cs_client */;
/*!50001 SET character_set_results     = @saved_cs_results */;
/*!50001 SET collation_connection      = @saved_col_connection */;

--
-- Final view structure for view `V_PLATE_RNAI`
--

/*!50001 DROP VIEW IF EXISTS `V_PLATE_RNAI`*/;
/*!50001 SET @saved_cs_client          = @@character_set_client */;
/*!50001 SET @saved_cs_results         = @@character_set_results */;
/*!50001 SET @saved_col_connection     = @@collation_connection */;
/*!50001 SET character_set_client      = utf8mb4 */;
/*!50001 SET character_set_results     = utf8mb4 */;
/*!50001 SET collation_connection      = utf8mb4_0900_ai_ci */;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`root`@`localhost` SQL SECURITY DEFINER */
/*!50001 VIEW `V_PLATE_RNAI` (`PLATE_ID`,`WELL`,`VENDOR_REAGENT_ID`,`CONC`) AS select `p`.`plate_id` AS `plate_id`,`c`.`well` AS `well`,`c`.`notebook_ref` AS `notebook_ref`,(`c`.`CONC` * 1e-3) AS `conc*1e-3` from ((`config` `c` join `plate` `p`) join `plating_sequence` `ps`) where ((`p`.`config_id` = `c`.`config_id`) and (`p`.`TYPE_ID` = `ps`.`TYPE_ID`) and (convert(`c`.`well` using utf8mb4) = convert(`ps`.`well` using utf8mb4)) and `c`.`config_id` in (select `c`.`config_id` from `config` `c` where `c`.`notebook_ref` in (select `r`.`vendor_reagent_id` from `sfl_assay`.`rnai` `r`))) order by `ps`.`SEQ` */;
/*!50001 SET character_set_client      = @saved_cs_client */;
/*!50001 SET character_set_results     = @saved_cs_results */;
/*!50001 SET collation_connection      = @saved_col_connection */;

--
-- Final view structure for view `v_plate`
--

/*!50001 DROP VIEW IF EXISTS `v_plate`*/;
/*!50001 SET @saved_cs_client          = @@character_set_client */;
/*!50001 SET @saved_cs_results         = @@character_set_results */;
/*!50001 SET @saved_col_connection     = @@collation_connection */;
/*!50001 SET character_set_client      = utf8mb4 */;
/*!50001 SET character_set_results     = utf8mb4 */;
/*!50001 SET collation_connection      = utf8mb4_0900_ai_ci */;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`root`@`localhost` SQL SECURITY DEFINER */
/*!50001 VIEW `v_plate` (`PLATE_ID`,`WELL`,`COMPOUND_ID`,`NOTEBOOK_REF`,`FORM`,`CONC`) AS select `p`.`plate_id` AS `plate_id`,`c`.`well` AS `well`,`c`.`compound_id` AS `compound_id`,`c`.`notebook_ref` AS `notebook_ref`,`c`.`FORM` AS `form`,`c`.`CONC` AS `conc` from ((`config` `c` join `plate` `p`) join `plating_sequence` `ps`) where ((`p`.`config_id` = `c`.`config_id`) and (`p`.`TYPE_ID` = `ps`.`TYPE_ID`) and (convert(`c`.`well` using utf8mb4) = convert(`ps`.`well` using utf8mb4))) order by `ps`.`SEQ` */;
/*!50001 SET character_set_client      = @saved_cs_client */;
/*!50001 SET character_set_results     = @saved_cs_results */;
/*!50001 SET collation_connection      = @saved_col_connection */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2022-03-20  8:02:53
