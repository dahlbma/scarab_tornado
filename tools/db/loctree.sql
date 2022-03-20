-- MySQL dump 10.13  Distrib 8.0.26, for Linux (x86_64)
--
-- Host: localhost    Database: loctree
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
-- Table structure for table `location_access`
--

DROP TABLE IF EXISTS `location_access`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `location_access` (
  `LOC_ID` text,
  `SEE` text,
  `USE` text,
  `CHANGE` text,
  `REMOVE` text
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `location_hierarchy`
--

DROP TABLE IF EXISTS `location_hierarchy`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `location_hierarchy` (
  `pkey` bigint NOT NULL AUTO_INCREMENT,
  `PARENT_TYPE` double DEFAULT NULL,
  `CHILD_TYPE` bigint DEFAULT NULL,
  PRIMARY KEY (`pkey`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `location_id_sequence`
--

DROP TABLE IF EXISTS `location_id_sequence`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `location_id_sequence` (
  `pkey` int NOT NULL,
  PRIMARY KEY (`pkey`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `location_type`
--

DROP TABLE IF EXISTS `location_type`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `location_type` (
  `type_id` int NOT NULL AUTO_INCREMENT,
  `NAME` text,
  `USE_SUBPOS` bigint DEFAULT NULL,
  `TERMINATED_DATE` text,
  `LABEL_FORMAT` text,
  `SUBPOS` double DEFAULT NULL,
  `APPLICATION_TAG` text,
  `PIXMAP` text,
  `PIXMAP_OPEN` text,
  PRIMARY KEY (`type_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `locations`
--

DROP TABLE IF EXISTS `locations`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `locations` (
  `loc_id` varchar(40) NOT NULL,
  `parent` varchar(40) DEFAULT NULL,
  `TYPE_ID` bigint DEFAULT NULL,
  `CREATED_DATE` datetime DEFAULT NULL,
  `CREATED_BY` bigint DEFAULT NULL,
  `NAME` text,
  `TERMINATED_DATE` text,
  `PRIORITY` double DEFAULT NULL,
  PRIMARY KEY (`loc_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `subpositions`
--

DROP TABLE IF EXISTS `subpositions`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `subpositions` (
  `TYPE_ID` bigint DEFAULT NULL,
  `SEQ` bigint DEFAULT NULL,
  `POS` text,
  UNIQUE KEY `subpositions_uk_type_seq` (`TYPE_ID`,`SEQ`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Temporary view structure for view `v_all_locations`
--

DROP TABLE IF EXISTS `v_all_locations`;
/*!50001 DROP VIEW IF EXISTS `v_all_locations`*/;
SET @saved_cs_client     = @@character_set_client;
/*!50503 SET character_set_client = utf8mb4 */;
/*!50001 CREATE VIEW `v_all_locations` AS SELECT 
 1 AS `LOC_ID`,
 1 AS `PARENT`,
 1 AS `TYPE_ID`,
 1 AS `NAME`,
 1 AS `path`,
 1 AS `PRIORITY`,
 1 AS `TERMINATED_DATE`,
 1 AS `level`*/;
SET character_set_client = @saved_cs_client;

--
-- Final view structure for view `v_all_locations`
--

/*!50001 DROP VIEW IF EXISTS `v_all_locations`*/;
/*!50001 SET @saved_cs_client          = @@character_set_client */;
/*!50001 SET @saved_cs_results         = @@character_set_results */;
/*!50001 SET @saved_col_connection     = @@collation_connection */;
/*!50001 SET character_set_client      = utf8mb4 */;
/*!50001 SET character_set_results     = utf8mb4 */;
/*!50001 SET collation_connection      = utf8mb4_0900_ai_ci */;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`root`@`localhost` SQL SECURITY DEFINER */
/*!50001 VIEW `v_all_locations` (`LOC_ID`,`PARENT`,`TYPE_ID`,`NAME`,`path`,`PRIORITY`,`TERMINATED_DATE`,`level`) AS with recursive `all_locs` as (select `locations`.`loc_id` AS `loc_id`,`locations`.`parent` AS `parent`,`locations`.`TYPE_ID` AS `type_id`,`locations`.`NAME` AS `name`,`locations`.`NAME` AS `path`,`locations`.`PRIORITY` AS `priority`,`locations`.`TERMINATED_DATE` AS `terminated_date`,1 AS `level` from `locations` where (`locations`.`parent` is null) union all select `this`.`loc_id` AS `loc_id`,`this`.`parent` AS `parent`,`this`.`TYPE_ID` AS `type_id`,`this`.`NAME` AS `name`,concat(`prior`.`path`,'/',`this`.`NAME`) AS `path`,`this`.`PRIORITY` AS `priority`,`this`.`TERMINATED_DATE` AS `terminated_date`,(`prior`.`level` + 1) AS `prior.level + 1` from (`all_locs` `prior` join `locations` `this` on((`this`.`parent` = `prior`.`loc_id`)))) select `e`.`loc_id` AS `loc_id`,`e`.`parent` AS `parent`,`e`.`type_id` AS `type_id`,`e`.`name` AS `name`,`e`.`path` AS `path`,`e`.`priority` AS `priority`,`e`.`terminated_date` AS `terminated_date`,`e`.`level` AS `level` from `all_locs` `e` order by `e`.`level` */;
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

-- Dump completed on 2022-03-20 13:20:33
