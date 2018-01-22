
# Author: Alfredo Sanchez-Tojar, MPIO (Seewiesen) and ICL (Silwood Park), alfredo.tojar@gmail.com
# Github profile: https://github.com/ASanchez-Tojar

# Script created on the 6th of July, 2016

########################################################################################################
# Description of script and Instructions
########################################################################################################

# This script is to prepare the databases needed to run heritability analyses
# on dominance on Lundy Island


########################################################################################################
# Packages needed
########################################################################################################

# packages needed to be loaded for this script

# libraries needed

library(pedantics)


# Clear memory and get to know where you are
rm(list=ls())
#getwd()


########################################################################################################
# Databases needed
########################################################################################################

# loading dominance data

#dominance<-read.table("001/rank.TLandM.VB.fitness.sex.csv",sep=",",header=TRUE)
dominance<-read.table("001/rank.TLandM.VB.fitness_sim.csv",sep=",",header=TRUE)
#Data <- as.data.frame(Data

# loading the Lundy pedigree

Ped <- read.table("Pedigree/PedigreeUpToIncl2016-mysterydissapear_versionwithNA.txt",sep="\t",header=TRUE)


########################################################################################################
# Creating variables and databases
########################################################################################################

# reducing dominance database

dominance <- dominance[,c("BirdID_eventSW","BirdID","StElo",
                          "eventSW","sex","cohort","age",
                          "tarsus","season","elo.z.event")]


# animal is needed for calling pedigree from MCMCglmm

dominance$animal <- dominance$BirdID


# reducing Ped to ID, mum and dad, easier to handle

Ped <- Ped[,c(1,2,3)]


# fixPedigree() from pedantics formats the pedigree in the way MCMCglmm wants it

Ped<-fixPedigree(Ped)


# reducing dominance database to those that we have pedigree data from

Ped.id <- Ped$id

dominance$BirdID <- as.numeric(dominance$BirdID)

dominance$inped <- ifelse((dominance$BirdID %in% Ped.id),
                          1,
                          0)

dominance <- subset(dominance,dominance$inped==1)


# Formatting some variables

dominance$animal <- as.factor(dominance$animal)
dominance$BirdID <- as.factor(dominance$BirdID)


# making id, dam and sire factors

for (x in 1:3) Ped[, x] <- as.factor(Ped[, x])



########################################################################################################
# Adding social Brood to each BirdID
########################################################################################################

# First, getting original BroodRef

# SELECT tblBirdID.BirdID, tblBirdID.BroodRef
# FROM tblBirdID;

BroodRef <- read.table("001/BirdID_BroodRef.csv",
                       sep=",",header=TRUE)

# removing the NA's from BroodRef as they don't help

BroodRef.2 <- BroodRef[!(is.na(BroodRef$BroodRef)),]


# Second, getting FosterBroodRef

# SELECT tblBirdID.BirdID, tblFosterBroods.FosterBrood
# FROM tblBirdID INNER JOIN tblFosterBroods ON tblBirdID.BirdID = tblFosterBroods.BirdID;

FosterBroodRef <- read.table("001/BirdID_FosterBroodRef.csv",
                             sep=",",header=TRUE)

# removing the NA's from FosterBroodRef as they don't have any information

FosterBroodRef.2 <- FosterBroodRef[!(is.na(FosterBroodRef$FosterBrood)) &
                                     !(is.na(FosterBroodRef$BirdID)),
                                   c("BirdID","FosterBrood")]



# adding BroodRef to dominance

dominance.brood <- merge(dominance,BroodRef.2,
                         by="BirdID",
                         all.x=TRUE)


# adding FosterBroodRef to dominance

dominance.brood.2 <- merge(dominance.brood,FosterBroodRef.2,
                           by="BirdID",
                           all.x=TRUE)


# creating a variable that says whether bird was crossfostered or not

dominance.brood.2$xfoster <- ifelse(!(is.na(dominance.brood.2$FosterBrood)),
                                    1,
                                    0)


# Generating the BroodRef where bird stayed from day 2 (after crossfoster) on

dominance.brood.2$SocialBroodRef<-dominance.brood.2$FosterBrood

dominance.brood.2$SocialBroodRef<-ifelse(is.na(dominance.brood.2$SocialBroodRef),
                                         dominance.brood.2$BroodRef,
                                         dominance.brood.2$SocialBroodRef)


########################################################################################################
# Adding social parents according to SocialBroodRef
########################################################################################################

# First, database until 2013.

# SELECT tblBroods.BroodRef, tblBroods.SocialDadID, tblBroods.SocialDadCertain, tblBroods.SocialMumID, tblBroods.SocialMumCertain
# FROM tblBroods LEFT JOIN tblBroodEvents ON tblBroods.BroodRef = tblBroodEvents.BroodRef
# WHERE (((tblBroodEvents.EventNumber)=0) AND ((tblBroodEvents.EventDate)<#1/1/2014#));

SocialParents <- read.table("001/BroodRef_SocialParents_until2013.csv",
                            sep=",",header=TRUE)


# Second I need to assemble the database from 2014 on


SocialDads <- read.table("001/BroodRef_SocialDadID_IDCertain_2014-2016_Updated.csv",
                         sep=",",header=TRUE)

SocialMums <- read.table("001/BroodRef_SocialMumID_IDCertain_2014-2016_Updated.csv",
                         sep=",",header=TRUE)

# There are no differences in BroodRef between both database

# setdiff(SocialDads$BroodRef,SocialMums$BroodRef)

# Therefore, I can merge them without losing information

SocialParents.2 <- merge(SocialDads,SocialMums,
                         by="BroodRef")



# Now, I can rbind() to get the final database to be added to dominance.brood.2


SocialParents.full <- rbind(SocialParents,SocialParents.2)


# assigning SocialCertain TRUE to those with missing SocialParents. Easier
# to have an idea what I miss than how it is coded now

SocialParents.full$SocialDadCertain <- ifelse(is.na(SocialParents.full$SocialDadID),
                                              TRUE,
                                              SocialParents.full$SocialDadCertain)

SocialParents.full$SocialMumCertain <- ifelse(is.na(SocialParents.full$SocialMumID),
                                              TRUE,
                                              SocialParents.full$SocialMumCertain)



# FINALLY adding SocialParents to each BirdID

dominance.brood.parents <- merge(dominance.brood.2,
                                 SocialParents.full,
                                 by.x="SocialBroodRef",
                                 by.y="BroodRef",
                                 all.x=TRUE)



########################################################################################################
# Adding NestboxName
########################################################################################################

# SELECT tblBroods.BroodRef, tblNestboxes.NestboxName
# FROM tblNestboxes RIGHT JOIN tblBroods ON tblNestboxes.NestboxRef = tblBroods.NestboxRef;

NestboxName <- read.table("001/BroodRef_NestboxName.csv",
                          sep=",",header=TRUE)


# now merging

dominance.brood.parents.nb <- merge(dominance.brood.parents,
                                    NestboxName,
                                    by.x="SocialBroodRef",
                                    by.y="BroodRef",
                                    all.x=TRUE)



########################################################################################################
# Adding sires and dams
########################################################################################################


h2dominance <- merge(dominance.brood.parents.nb,
                     Ped,
                     by.x="BirdID",
                     by.y="id",
                     all.x=TRUE)


h2dominance <- h2dominance[,c("BirdID","sex","cohort","tarsus",
                              "BirdID_eventSW","eventSW","season","StElo",
                              "elo.z.event","age",
                              "BroodRef","FosterBrood","SocialBroodRef","xfoster",
                              "SocialDadID","SocialDadCertain",
                              "SocialMumID","SocialMumCertain",
                              "sire","dam","animal","inped",
                              "NestboxName")]


write.csv(h2dominance,"002/h2dominance.csv",row.names=FALSE)



# # loading the database that contains social parents and social brood per BirdID
# 
# family <- read.table("BirdID_Parents_and_broodID_20160705.csv",sep=",",header=TRUE)
# 
# 
# # merging Data with family to put all individual information together
# 
# Data.v2 <- merge(Data, family, by="BirdID",all.x=TRUE)
# 
# 
# # Now I want to add the genetic mother and father for each BirdID. I can use the pedigree for this.
# 
# Data.v3 <- merge(Data.v2, Ped, by.x="BirdID",by.y="id",all.x=TRUE)
# 
# 
# # adding the nestbox to the database
# 
# nestbox <- read.table("brood_and_nestbox.csv",sep=",",header=TRUE)
# 
# 
# Data.v3 <- merge(Data.v3,nestbox,by="BroodName",all.x=TRUE)
# 
# 
# # removing NA's for the variables to be used in this model
# 
# Data2.1 <- subset(Data, !(is.na(Data$SexEstimate)))
# Data2.1 <- subset(Data2.1, !(is.na(Data2.1$Tarsus.mean)))
# 
# # removing NA's from the fixed effects and broodID
# 
# Data.v4 <- subset(Data.v3, !(is.na(Data.v3$SexEstimate)))
# Data.v5 <- subset(Data.v4, !(is.na(Data.v4$Tarsus.mean)))
# Data.v6 <- subset(Data.v5, !(is.na(Data.v5$BroodName)))
# Data.v7 <- subset(Data.v6,Data.v6$BroodName!="")
# 
# # removing NA's from SocialMumID and removing IDcertain=no
# 
# Data.v8 <- subset(Data.v5, !(is.na(Data.v5$SocialMumID)))
# Data.v9 <- subset(Data.v8,Data.v8$SocialMumCertain!=FALSE)
# 
# # removing NA's from SocialMumID and removing IDcertain=no
# 
# Data.v10 <- subset(Data.v5, !(is.na(Data.v5$SocialDadID)))
# Data.v11 <- subset(Data.v10,Data.v10$SocialDadCertain!=FALSE)
# 
# # removing NA's from SocialMumID and removing IDcertain=no
# 
# Data.v12 <- subset(Data.v5, !(is.na(Data.v5$dam)))
# 
# 
# # removing
# 
# Data.v13 <- subset(Data.v7, !(is.na(Data.v7$NestboxName)))

