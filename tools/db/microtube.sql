-- MySQL dump 10.13  Distrib 8.0.26, for Linux (x86_64)
--
-- Host: localhost    Database: microtube
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
-- Table structure for table `matrix`
--

DROP TABLE IF EXISTS `matrix`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `matrix` (
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `matrix_id` varchar(40) DEFAULT NULL,
  `CREATED_DATE` datetime DEFAULT NULL,
  `location` varchar(80) DEFAULT NULL,
  `COMMENTS` text,
  PRIMARY KEY (`pkey`),
  UNIQUE KEY `matrix_id_ix` (`matrix_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `matrix_position`
--

DROP TABLE IF EXISTS `matrix_position`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `matrix_position` (
  `seq` int NOT NULL,
  `position` varchar(10) DEFAULT NULL,
  PRIMARY KEY (`seq`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `matrix_tube`
--

DROP TABLE IF EXISTS `matrix_tube`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `matrix_tube` (
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `matrix_id` varchar(40) DEFAULT NULL,
  `position` varchar(10) DEFAULT NULL,
  `tube_id` varchar(10) DEFAULT NULL,
  PRIMARY KEY (`pkey`),
  KEY `tube_id` (`tube_id`),
  KEY `matrix_id` (`matrix_id`),
  CONSTRAINT `matrix_tube_ibfk_1` FOREIGN KEY (`tube_id`) REFERENCES `tube` (`tube_id`),
  CONSTRAINT `matrix_tube_ibfk_2` FOREIGN KEY (`matrix_id`) REFERENCES `matrix` (`matrix_id`)
) ENGINE=InnoDB AUTO_INCREMENT=127 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `matrix_tube_invalid`
--

DROP TABLE IF EXISTS `matrix_tube_invalid`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `matrix_tube_invalid` (
  `PKEY` bigint DEFAULT NULL,
  `matrix_id` varchar(40) DEFAULT NULL,
  `position` varchar(10) DEFAULT NULL,
  `tube_id` varchar(10) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tube`
--

DROP TABLE IF EXISTS `tube`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `tube` (
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `tube_id` varchar(10) DEFAULT NULL,
  `notebook_ref` varchar(16) DEFAULT NULL,
  `VOLUME` double DEFAULT NULL,
  `CONC` double DEFAULT NULL,
  `TDATE` datetime DEFAULT NULL,
  `SOURCE` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `CREATED_BY_PKEY` bigint DEFAULT NULL,
  `UPDATED_DATE` datetime DEFAULT NULL,
  `UPDATED_BY_PKEY` double DEFAULT NULL,
  `location` varchar(60) DEFAULT NULL,
  `THAW_COUNT` double DEFAULT NULL,
  PRIMARY KEY (`pkey`),
  UNIQUE KEY `tube_id_ix` (`tube_id`),
  KEY `notebook_ref` (`notebook_ref`),
  CONSTRAINT `tube_ibfk_1` FOREIGN KEY (`notebook_ref`) REFERENCES `bcpvs`.`batch` (`notebook_ref`)
) ENGINE=InnoDB AUTO_INCREMENT=40 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tube_changes`
--

DROP TABLE IF EXISTS `tube_changes`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `tube_changes` (
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `tube_id` varchar(10) DEFAULT NULL,
  `notebook_ref` varchar(30) DEFAULT NULL,
  `VOLUME` text,
  `DELTAVOLUME` text,
  `CONC` text,
  `TDATE` text,
  `SOURCE` text,
  `CREATED_DATE` text,
  `CREATED_BY_PKEY` text,
  `location` varchar(60) DEFAULT NULL,
  PRIMARY KEY (`pkey`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tube_invalid`
--

DROP TABLE IF EXISTS `tube_invalid`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `tube_invalid` (
  `PKEY` bigint DEFAULT NULL,
  `tube_id` varchar(10) DEFAULT NULL,
  `notebook_ref` varchar(16) DEFAULT NULL,
  `VOLUME` double DEFAULT NULL,
  `CONC` double DEFAULT NULL,
  `TDATE` datetime DEFAULT NULL,
  `SOURCE` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `CREATED_BY_PKEY` bigint DEFAULT NULL,
  `UPDATED_DATE` datetime DEFAULT NULL,
  `UPDATED_BY_PKEY` double DEFAULT NULL,
  `location` varchar(60) DEFAULT NULL,
  `THAW_COUNT` double DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Temporary view structure for view `v_matrix`
--

DROP TABLE IF EXISTS `v_matrix`;
/*!50001 DROP VIEW IF EXISTS `v_matrix`*/;
SET @saved_cs_client     = @@character_set_client;
/*!50503 SET character_set_client = utf8mb4 */;
/*!50001 CREATE VIEW `v_matrix` AS SELECT 
 1 AS `MATRIX_ID`,
 1 AS `CREATED_DATE`,
 1 AS `LOC_ID`,
 1 AS `LOCATION`,
 1 AS `COMMENTS`*/;
SET character_set_client = @saved_cs_client;

--
-- Temporary view structure for view `v_matrix_locations`
--

DROP TABLE IF EXISTS `v_matrix_locations`;
/*!50001 DROP VIEW IF EXISTS `v_matrix_locations`*/;
SET @saved_cs_client     = @@character_set_client;
/*!50503 SET character_set_client = utf8mb4 */;
/*!50001 CREATE VIEW `v_matrix_locations` AS SELECT 
 1 AS `LOC_ID`,
 1 AS `PATH`*/;
SET character_set_client = @saved_cs_client;

--
-- Temporary view structure for view `v_matrix_tube`
--

DROP TABLE IF EXISTS `v_matrix_tube`;
/*!50001 DROP VIEW IF EXISTS `v_matrix_tube`*/;
SET @saved_cs_client     = @@character_set_client;
/*!50503 SET character_set_client = utf8mb4 */;
/*!50001 CREATE VIEW `v_matrix_tube` AS SELECT 
 1 AS `MATRIX_ID`,
 1 AS `POSITION`,
 1 AS `TUBE_ID`*/;
SET character_set_client = @saved_cs_client;

--
-- Final view structure for view `v_matrix`
--

/*!50001 DROP VIEW IF EXISTS `v_matrix`*/;
/*!50001 SET @saved_cs_client          = @@character_set_client */;
/*!50001 SET @saved_cs_results         = @@character_set_results */;
/*!50001 SET @saved_col_connection     = @@collation_connection */;
/*!50001 SET character_set_client      = utf8mb4 */;
/*!50001 SET character_set_results     = utf8mb4 */;
/*!50001 SET collation_connection      = utf8mb4_0900_ai_ci */;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`root`@`localhost` SQL SECURITY DEFINER */
/*!50001 VIEW `v_matrix` (`MATRIX_ID`,`CREATED_DATE`,`LOC_ID`,`LOCATION`,`COMMENTS`) AS select `m`.`matrix_id` AS `MATRIX_ID`,`m`.`CREATED_DATE` AS `CREATED_DATE`,`m`.`location` AS `LOCATION`,ifnull(convert(`l`.`PATH` using utf8mb4),convert(`m`.`location` using utf8mb4)) AS `location`,`m`.`COMMENTS` AS `COMMENTS` from (`matrix` `m` left join `v_matrix_locations` `l` on((convert(`m`.`location` using utf8mb4) = convert(`l`.`LOC_ID` using utf8mb4)))) */;
/*!50001 SET character_set_client      = @saved_cs_client */;
/*!50001 SET character_set_results     = @saved_cs_results */;
/*!50001 SET collation_connection      = @saved_col_connection */;

--
-- Final view structure for view `v_matrix_locations`
--

/*!50001 DROP VIEW IF EXISTS `v_matrix_locations`*/;
/*!50001 SET @saved_cs_client          = @@character_set_client */;
/*!50001 SET @saved_cs_results         = @@character_set_results */;
/*!50001 SET @saved_col_connection     = @@collation_connection */;
/*!50001 SET character_set_client      = utf8mb4 */;
/*!50001 SET character_set_results     = utf8mb4 */;
/*!50001 SET collation_connection      = utf8mb4_0900_ai_ci */;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`root`@`localhost` SQL SECURITY DEFINER */
/*!50001 VIEW `v_matrix_locations` (`LOC_ID`,`PATH`) AS select `loctree`.`v_all_locations`.`LOC_ID` AS `loc_id`,`loctree`.`v_all_locations`.`path` AS `PATH` from `loctree`.`v_all_locations` where `loctree`.`v_all_locations`.`TYPE_ID` in (select `loctree`.`location_hierarchy`.`PARENT_TYPE` from `loctree`.`location_hierarchy` where (`loctree`.`location_hierarchy`.`CHILD_TYPE` = 20)) */;
/*!50001 SET character_set_client      = @saved_cs_client */;
/*!50001 SET character_set_results     = @saved_cs_results */;
/*!50001 SET collation_connection      = @saved_col_connection */;

--
-- Final view structure for view `v_matrix_tube`
--

/*!50001 DROP VIEW IF EXISTS `v_matrix_tube`*/;
/*!50001 SET @saved_cs_client          = @@character_set_client */;
/*!50001 SET @saved_cs_results         = @@character_set_results */;
/*!50001 SET @saved_col_connection     = @@collation_connection */;
/*!50001 SET character_set_client      = utf8mb4 */;
/*!50001 SET character_set_results     = utf8mb4 */;
/*!50001 SET collation_connection      = utf8mb4_0900_ai_ci */;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`root`@`localhost` SQL SECURITY DEFINER */
/*!50001 VIEW `v_matrix_tube` (`MATRIX_ID`,`POSITION`,`TUBE_ID`) AS select `m`.`matrix_id` AS `MATRIX_ID`,`m`.`position` AS `POSITION`,`m`.`tube_id` AS `TUBE_ID` from (`matrix_tube` `m` join `matrix_position` `p`) where (convert(`m`.`position` using utf8mb4) = convert(`p`.`position` using utf8mb4)) order by `m`.`matrix_id`,`p`.`seq` */;
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

-- Dump completed on 2022-03-20 13:20:51
