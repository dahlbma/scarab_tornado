-- MySQL dump 10.13  Distrib 8.0.26, for Linux (x86_64)
--
-- Host: localhost    Database: glass
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
-- Table structure for table `box_sequence`
--

DROP TABLE IF EXISTS `box_sequence`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `box_sequence` (
  `coordinate` varchar(8) COLLATE utf8_unicode_ci DEFAULT NULL,
  UNIQUE KEY `box_pos_idx` (`coordinate`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `globals`
--

DROP TABLE IF EXISTS `globals`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `globals` (
  `pkey` bigint NOT NULL,
  `VARIABLE_NAME` text,
  `VALUE` text,
  PRIMARY KEY (`pkey`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `hive_stats`
--

DROP TABLE IF EXISTS `hive_stats`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `hive_stats` (
  `pkey` bigint NOT NULL,
  `user_id` varchar(20) DEFAULT NULL,
  `TYPE` text,
  `TDATE` datetime DEFAULT NULL,
  `HOSTNAME` text,
  PRIMARY KEY (`pkey`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `vial`
--

DROP TABLE IF EXISTS `vial`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `vial` (
  `vial_id` varchar(7) NOT NULL,
  `TYPE_ID` bigint DEFAULT NULL,
  `notebook_ref` varchar(16) DEFAULT NULL,
  `FORM` text,
  `CONC` double DEFAULT NULL,
  `location` varchar(40) DEFAULT NULL,
  `COMMENTS` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `CREATOR_PKEY` bigint DEFAULT NULL,
  `UPDATED_DATE` datetime DEFAULT NULL,
  `UPDATER_PKEY` double DEFAULT NULL,
  `TARE` double DEFAULT NULL,
  `GROSS` double DEFAULT NULL,
  `NET` double DEFAULT NULL,
  `pos` varchar(8) DEFAULT NULL,
  PRIMARY KEY (`vial_id`),
  UNIQUE KEY `vial_idx` (`vial_id`),
  KEY `notebook_ref` (`notebook_ref`),
  KEY `loc_idx` (`location`),
  KEY `pos_idx` (`pos`),
  CONSTRAINT `vial_ibfk_1` FOREIGN KEY (`notebook_ref`) REFERENCES `bcpvs`.`batch` (`notebook_ref`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `vial_invalid`
--

DROP TABLE IF EXISTS `vial_invalid`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `vial_invalid` (
  `vial_id` varchar(7) DEFAULT NULL,
  `TYPE_ID` bigint DEFAULT NULL,
  `NOTEBOOK_REF` text,
  `FORM` text,
  `CONC` double DEFAULT NULL,
  `location` varchar(40) DEFAULT NULL,
  `COMMENTS` text,
  `CREATED_DATE` datetime DEFAULT NULL,
  `CREATOR_PKEY` bigint DEFAULT NULL,
  `UPDATED_DATE` datetime DEFAULT NULL,
  `UPDATER_PKEY` double DEFAULT NULL,
  `TARE` double DEFAULT NULL,
  `GROSS` double DEFAULT NULL,
  `NET` double DEFAULT NULL,
  `pos` varchar(8) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `vial_log`
--

DROP TABLE IF EXISTS `vial_log`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `vial_log` (
  `vial_id` varchar(7) DEFAULT NULL,
  `UPDATER_PKEY` bigint DEFAULT NULL,
  `UPDATED_DATE` datetime DEFAULT NULL,
  `CHANGES` text,
  KEY `vial_id` (`vial_id`),
  CONSTRAINT `vial_log_ibfk_1` FOREIGN KEY (`vial_id`) REFERENCES `vial` (`vial_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `vial_log_invalid`
--

DROP TABLE IF EXISTS `vial_log_invalid`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `vial_log_invalid` (
  `vial_id` varchar(7) DEFAULT NULL,
  `UPDATER_PKEY` bigint DEFAULT NULL,
  `UPDATED_DATE` datetime DEFAULT NULL,
  `CHANGES` text
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

-- Dump completed on 2022-03-20 13:20:19
