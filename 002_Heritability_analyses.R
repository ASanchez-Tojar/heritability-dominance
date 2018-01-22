
# Author: Alfredo Sanchez-Tojar, MPIO (Seewiesen) and ICL (Silwood Park), alfredo.tojar@gmail.com
# Github profile: https://github.com/ASanchez-Tojar

# Script created on the 6th of July, 2016

########################################################################################################
# Description of script and Instructions
########################################################################################################

# This script is analyze dominance heritability on Lundy Island

# Following MCMCglmm tutorial from http://www.wildanimalmodels.org/


########################################################################################################
# Packages needed
########################################################################################################

# packages needed to be loaded for this script

# libraries needed

library(MCMCglmm)
library(pedantics)
library(forestplot)
library(plyr)


# Clear memory and get to know where you are
rm(list=ls())
#getwd()


########################################################################################################
# Databases needed
########################################################################################################

# loading h2dominance data

h2dominance <- read.table("002/h2dominance.csv",sep=",",header=TRUE)


# loading the Lundy pedigree

Ped <- read.table("Pedigree/PedigreeUpToIncl2016-mysterydissapear_versionwithNA.txt",sep="\t",header=TRUE)

Ped.fix<-fixPedigree(Ped)

dataavailable <- ifelse(Ped.fix$id %in% h2dominance$BirdID,1,0)

Ped.fix.2 <- cbind(Ped.fix,dataavailable)

BirdIDcohort <- read.table("Pedigree/BirdID_Cohort_20180106.csv",
                           header=TRUE,sep=",")

names(BirdIDcohort) <- c("id","Cohort")

Ped.fix.3 <- join(Ped.fix.2, BirdIDcohort,by="id",type = "left")

#Removing two birds with no cohort, and no parents
#Ped.fix.3[is.na(Ped.fix.3$Cohort),]
Ped.fix.4 <- Ped.fix.3[!(is.na(Ped.fix.3$Cohort)),]

prunedPed <- prunePed(Ped.fix.4,h2dominance$BirdID)


# tiff("Pedigree.tiff", height=40, width=50,
#      units='cm', compression="lzw", res=300)     
# 
# rgbing <- c(255,255,255)
# darkblue <- c(31,120,180)/rgbing
# chocolate1 <- c(255,127,36)/rgbing
# 
# sexcolours<-c(rgb(chocolate1[1],chocolate1[2],chocolate1[3],0.5),
#               rgb(darkblue[1],darkblue[2],darkblue[3],0.5))
# 
# drawPedigree(Ped.fix.4[,c("id","dam","sire")],
#              cohorts = Ped.fix.4$Cohort,
#              writeCohortLabels='y',
#              cohortLabs.cex = 1.5,
#              sexColours = sexcolours,
#              #retain="pruned",
#              dataDots='y', 
#              dotSize = 0.003,
#              dat=Ped.fix.4$dataavailable)
#              #dots='y')
# 
# dev.off()

# # Pedigree stats
# stats.g<-pedigreeStats(Ped.fix.4[,c("id","dam","sire")],
#                        graphicalReport='n')
# 
# stats.g.pruned<-pedigreeStats(Ped.fix.4[,c("id","dam","sire")],
#                               dat=Ped.fix.4$dataavailable,
#                               graphicalReport='n')
# 
# stats.g.pruned.2<-pedigreeStats(prunedPed[,c("id","dam","sire")],
#                                 graphicalReport='n')
# 
# round(pedStatSummary(stats.g.pruned.2),4)
# round(pedStatSummary(stats.g,stats.g.pruned.2),4)
# #
# sink("full_and_pruned_Pedigree_stats.txt")
# round(pedStatSummary(stats.g,stats.g.pruned.2),4)
# sink()


# ######################################################################################################
# # MODEL 0: StElo ~ 1, random = ~ BirdID
# ######################################################################################################
# 
# # defining priors
# 
# prior0.1<-list(R = list(V = diag(1), nu = 0.002), 
#                G =list(G1 = list(V = diag(1), nu = 0.002)))
# 
# 
# #  model
# 
# model0.1 <- MCMCglmm(StElo ~ 1,
#                      random = ~BirdID, 
#                      pedigree = Ped,
#                      data = h2dominance, 
#                      prior = prior0.1,
#                      nitt = 7500000, 
#                      thin = 3600, 
#                      burnin = 750000,
#                      verbose=FALSE)
# 
# # save model if you haven't done it yet
# 
# save("model0.1", file = "mod0.1_StElo_BirdID_nodisplacements_7.5m_20160912.rda")
# 
# 
# # # loading the already run model
# # 
# # load("mod0_mod0.1_StElo_repeatability_aggression_4million_3600thin-20160705.rda")
# # summary(model0.1)
# 
# 
# # evaluating how confident we can be that MCMCglmm found good answers
# 
# # checking how the model ranplot(model0.1$Sol)
# 
# plot(model0.1$VCV)
# 
# 
# # checking if the chain was more or less stable in the value
# 
# #geweke.diag(model0.1$VCV, frac1=0.1, frac2=0.5)
# 
# geweke.plot(model0.1$VCV, frac1=0.1,
#             frac2=0.5, nbins = 20,
#             pvalue = 0.05, auto.layout = TRUE)
# 
# 
# #raftery.diag(model0.1$VCV, q=0.05, r=0.01, s=0.95, converge.eps=0.001)
# 
# 
# # checking if there is autocorrelation between the values predicted
# 
# autocorr(model0.1$VCV)
# 
# autocorr.plot(model0.1$VCV)
# 
# #heidel.diag(model0.1$Sol)
# #heidel.diag(model0.1$VCV)
# 
# 
# # checking the results
# 
# # variance explained by BirdID
# 
# posterior.mode(model0.1$VCV)
# 
# HPDinterval(model0.1$VCV)
# 
# 
# # estimating the percentage of variance explain by BirdID, i.e. checking Vi
# 
# model0.1.VP <- model0.1$VCV[,"BirdID"]+model0.1$VCV[,"units"]
# 
# posterior.mode(model0.1$VCV[,"BirdID"]/model0.1.VP)
# 
# 
# # plotting the results as a forest plot
# 
# # list of values to be plotted
# 
# Modoff <- structure(c(posterior.mode(model0.1$VCV[,"BirdID"]/model0.1.VP), # mode repeatability
#                       posterior.mode(model0.1$VCV[,"units"]/model0.1.VP), # mode residual
#                       HPDinterval(model0.1$VCV)[1,1]/posterior.mode(model0.1.VP), # lower CrI repeatability
#                       HPDinterval(model0.1$VCV)[2,1]/posterior.mode(model0.1.VP), # lower CrI residuals
#                       HPDinterval(model0.1$VCV)[1,2]/posterior.mode(model0.1.VP), # upper CrI repeatability
#                       HPDinterval(model0.1$VCV)[2,2]/posterior.mode(model0.1.VP)  # upper CrI residuals 
# ), .Dim=c(2L, 3L), .Dimnames = list(c("repeatability", "residuals"),
#                                     c("post.mode", "lower 95CrI", "upper 95CrI")))                    
# 
# 
# plot(1, 1,type="n",xlab="",ylab= "",xaxt="n",yaxt="n", bty="n")# using this nonsense because otherwise it overwrites existing plot
# 
# 
# 
# forestplot(c(" Vi \n----- \n Vp"," Vr \n----- \n Vp"),
#            mean = Modoff[, "post.mode"], 
#            lower = Modoff[, "lower 95CrI"], 
#            upper = Modoff[, "upper 95CrI"],
#            boxsize = 0.05,
#            xticks=c(0, 0.1, 0.2, 0.3, 0.4,0.5,0.6,0.7,0.8,0.9,1.0), line.margin = unit(50,"mm"), 
#            lineheight =  unit(6, "cm"),
#            grid = TRUE,
#            lwd.xaxis = 2,lwd.zero = 2, lwd.ci = 2,
#            col=fpColors(box = "black", lines = "black",
#                         zero = "lightgray", text = "black", axes = "black",
#                         hrz_lines = "black"),
#            title="Standardized Elo-rating based on aggressive interactions\nVi=birdID, Vr=residuals, Vp=Vi+Vr")


######################################################################################################
# MODEL 1: StElo ~ 1, random = ~ animal + BirdID
######################################################################################################

ptm <- proc.time()

# defining priors

prior1.1<-list(R = list(V = diag(1), nu = 0.002), 
               G =list(G1 = list(V = diag(1), nu = 0.002), 
                       G2 = list(V = diag(1), nu = 0.002)))


#h2dominance.1 <- h2dominance[h2dominance$inped==1,]


#  model

# model1.1 <- MCMCglmm(StElo ~ 1,
#                      random = ~animal+BirdID, 
#                      pedigree = prunedPed[,c(1,2,3)],
#                      data = h2dominance, 
#                      prior = prior1.1,
#                      nitt = 5000000,
#                      thin = 2500,
#                      burnin = 1250000,
#                      # nitt = 3000000,
#                      # thin = 3000,
#                      # burnin = 300000,
#                      verbose=FALSE)
# 
# proc.time() - ptm

# save model if you haven't done it yet

#save("model1.1", file = "models/mod1_h2_dominance.rda")
#save("model1.1", file = "models/mod1_h2_dominance_3m.rda")


# loading the already run model

load("models/mod1_h2_dominance.rda")
summary(model1.1)


# evaluating how confident we can be that MCMCglmm found good answers

# checking how the model ran

plot(model1.1$Sol)
plot(model1.1$VCV)


# checking if the chain was more or less stable in the value

#geweke.diag(model1.1$VCV, frac1=0.1, frac2=0.5)

geweke.plot(model1.1$VCV, frac1=0.1,
            frac2=0.5, nbins = 20,
            pvalue = 0.05, auto.layout = TRUE)


#raftery.diag(model1.1$VCV, q=0.05, r=0.01, s=0.95, converge.eps=0.001)


# checking if there is autocorrelation between the values predicted

autocorr(model1.1$Sol)

autocorr(model1.1$VCV)

autocorr.plot(model1.1$Sol)

autocorr.plot(model1.1$VCV)

#heidel.diag(model1.1$Sol)
#heidel.diag(model1.1$VCV)


# checking the results

# variance explained by the pedigree and BirdID

posterior.mode(model1.1$VCV)

HPDinterval(model1.1$VCV)


# estimating the percentage of variance explain by the pedigree and BirdID, i.e. checking Va and Vpe

model1.1.VP <- model1.1$VCV[,"animal"]+model1.1$VCV[,"BirdID"]+model1.1$VCV[,"units"]

model1.1.IDplusVA <- model1.1$VCV[,"animal"]+model1.1$VCV[,"BirdID"]

posterior.mode(model1.1.IDplusVA/model1.1.VP)

# raw variance
posterior.mode(model1.1$VCV[,"animal"])
round(mean(model1.1$VCV[,"animal"]),0)
HPDinterval(model1.1$VCV[,"animal"])

posterior.mode(model1.1$VCV[,"BirdID"])
round(mean(model1.1$VCV[,"BirdID"]),0)
HPDinterval(model1.1$VCV[,"BirdID"])

posterior.mode(model1.1$VCV[,"units"])
round(mean(model1.1$VCV[,"units"]),0)
HPDinterval(model1.1$VCV[,"units"])

# proportion of variance
posterior.mode(model1.1$VCV[,"animal"]/model1.1.VP)
round(mean(model1.1$VCV[,"animal"]/model1.1.VP),3)
HPDinterval(model1.1$VCV[,"animal"]/model1.1.VP)

posterior.mode(model1.1$VCV[,"BirdID"]/model1.1.VP)
round(mean(model1.1$VCV[,"BirdID"]/model1.1.VP),3)
HPDinterval(model1.1$VCV[,"BirdID"]/model1.1.VP)

posterior.mode(model1.1$VCV[,"units"]/model1.1.VP)
round(mean(model1.1$VCV[,"units"]/model1.1.VP),3)
HPDinterval(model1.1$VCV[,"units"]/model1.1.VP)



#########################################################
# plotting posterior distributions of heritability and PE
#########################################################

tiff("posterior_distributions.tiff", 
     height=15, width=30,
     units='cm', compression="lzw", res=600)

m <- rbind(c(1,2))


layout(m)

op <- par(oma = c(4,3,0.5,0.5) + 0.1,
          mar = c(1.5,3,1,0) + 0.1)

d.animal <- density(model1.1$VCV[,"animal"]/model1.1.VP)

plot(d.animal, 
     #type="n",
     xlab="",
     ylab= "",
     main="",
     xaxt="n",yaxt="n",xlim=c(0,0.5),ylim=c(0,6),
     family="serif")

axis(1,at=seq(0,0.5,0.1),
     las=1,
     cex.axis=1.5,
     family="serif") 

axis(2,at=seq(0,6,1),
     cex.axis=1.5,
     las=2,
     family="serif")

title(ylab = "density",
      line=0,
      outer = TRUE, cex.lab=3.25)

text(0.5,6,"(a)")

polygon(d.animal, col="grey75", border="black")
# hist(model1.1$VCV[,"animal"]/model1.1.VP,
#      add=TRUE,freq=FALSE,breaks=16)

lines(c(mean(model1.1$VCV[,"animal"]/model1.1.VP),
        mean(model1.1$VCV[,"animal"]/model1.1.VP)),
      c(0,6.2),lwd=2)

d.BirdID <- density(model1.1$VCV[,"BirdID"]/model1.1.VP)

plot(d.BirdID, 
     #type="n",
     main="",
     xlab="",
     ylab= "",
     xaxt="n",yaxt="n",xlim=c(0,0.5),ylim=c(0,6),
     family="serif")

axis(1,at=seq(0,0.5,0.1),
     las=1,
     cex.axis=1.5,
     family="serif") 

axis(2,at=seq(0,6,1),
     cex.axis=1.5,
     las=2,
     family="serif")

polygon(d.BirdID, col="grey75", border="black")

lines(c(mean(model1.1$VCV[,"BirdID"]/model1.1.VP),
        mean(model1.1$VCV[,"BirdID"]/model1.1.VP)),
      c(0,6.2),lwd=2)

title(xlab = "proportion of variance",
      outer = TRUE,
      line=2,
      cex.lab=3.25)

text(0.5,6,"(b)")

dev.off()

# hist(model1.1$VCV[,"BirdID"]/model1.1.VP,add=TRUE)


# plotting the results as a forest plot

# list of values to be plotted

# Modoff <- structure(c(posterior.mode(model1.1$VCV[,"animal"]/model1.1.VP), # mode heritability
#                       posterior.mode(model1.1.IDplusVA/model1.1.VP), # mode repeatability
#                       posterior.mode(model1.1$VCV[,"units"]/model1.1.VP), # mode residual
#                       HPDinterval(model1.1$VCV)[1,1]/posterior.mode(model1.1.VP), # lower CrI heritability
#                       (HPDinterval(model1.1$VCV)[1,1]+HPDinterval(model1.1$VCV)[2,1])/posterior.mode(model1.1.VP), # lower CrI repeatability
#                       HPDinterval(model1.1$VCV)[3,1]/posterior.mode(model1.1.VP), # lower CrI residuals
#                       HPDinterval(model1.1$VCV)[1,2]/posterior.mode(model1.1.VP), # upper CrI heritability
#                       (HPDinterval(model1.1$VCV)[1,2]+HPDinterval(model1.1$VCV)[2,2])/posterior.mode(model1.1.VP), # upper CrI repeatability
#                       HPDinterval(model1.1$VCV)[3,2]/posterior.mode(model1.1.VP)  # upper CrI residuals   
#                       
#                       
# ), .Dim=c(3L, 3L), .Dimnames = list(c("heritability", "repeatability", "residuals"),
#                                     c("post.mode", "lower 95CrI", "upper 95CrI")))
# 
# 
# Modoff <- structure(c(posterior.mode(model1.1$VCV[,"animal"]/model1.1.VP), # mode heritability
#                       posterior.mode(model1.1$VCV[,"BirdID"]/model1.1.VP), # mode repeatability
#                       HPDinterval(model1.1$VCV)[1,1]/posterior.mode(model1.1.VP), # lower CrI heritability
#                       HPDinterval(model1.1$VCV)[2,1]/posterior.mode(model1.1.VP),
#                       HPDinterval(model1.1$VCV)[1,2]/posterior.mode(model1.1.VP), # upper CrI heritability
#                       HPDinterval(model1.1$VCV)[2,2]/posterior.mode(model1.1.VP)
#                       
# ), .Dim=c(2L, 3L), .Dimnames = list(c("h2", "PE"),
#                                     c("post.mode", "lower 95CrI", "upper 95CrI")))
# 
# 
# plot(1, 1,type="n",xlab="",ylab= "",xaxt="n",yaxt="n", bty="n")# using this nonsense because otherwise it overwrites existing plot
# 
# 
# 
# forestplot(c("h2","PE"),
#            mean = Modoff[, "post.mode"], 
#            lower = Modoff[, "lower 95CrI"], 
#            upper = Modoff[, "upper 95CrI"],
#            boxsize = 0.1,
#            txt_gp = fpTxtGp(label = gpar(fontfamily = "serif",cex=2),
#                             ticks = gpar(fontfamily = "serif", cex=1.5),
#                             xlab  = gpar(fontfamily = "serif", cex = 2.5)),
#            xlab="proportion of variance",
#            xticks=c(0, 0.1, 0.2, 0.3, 0.4), line.margin = unit(50,"mm"), 
#            lineheight =  unit(8, "cm"),
#            grid = TRUE,
#            lwd.xaxis = 2,lwd.zero = 2, lwd.ci = 2,
#            col=fpColors(box = "dodgerblue4", lines = "dodgerblue4",
#                         zero = "black", text = "black", axes = "black",
#                         hrz_lines = "black"))


# ######################################################################################################
# # MODEL 2: StElo ~ sex + age + tarsus, random = ~ animal + BirdID
# ######################################################################################################
# 
# 
# # defining priors
# 
# prior2.1<-list(R = list(V = diag(1), nu = 0.002), 
#                G =list(G1 = list(V = diag(1), nu = 0.002), G2 = list(V = diag(1), nu = 0.002)))
# 
# 
# h2dominance.2 <- h2dominance[!(is.na(h2dominance$sex)) &
#                                !(is.na(h2dominance$age)) &
#                                h2dominance$inped==1,]
# 
# 
# #  model
# 
# model2.1 <- MCMCglmm(StElo ~ scale(sex)+
#                      scale(age),
#                      random = ~animal+BirdID, 
#                      pedigree = Ped,
#                      data = h2dominance.2, 
#                      prior = prior2.1,
#                      nitt = 5000000, 
#                      thin = 2500, 
#                      burnin = 1250000,
#                      verbose=FALSE)
# 
# # save model if you haven't done it yet
# #save("model1.1", file = "models/mod1_h2_dominance_fixedfactors.rda")
# 
# 
# # # loading the already run model
# # 
# # load("mod2.1_StElo_h2_fixedeffects_nodisplacments_7.5m_20160909.rda")
# # summary(model2.1)
# 
# load("mod2.1_StElo_h2_nodisplacments_7.5m_20160910.rda")
# 
# 
# 
# # evaluating how confident we can be that MCMCglmm found good answers
# 
# # checking how the model ran
# 
# plot(model2.1$Sol)
# plot(model2.1$VCV)
# 
# 
# # checking if the chain was more or less stable in the value
# 
# #geweke.diag(model2.1$Sol, frac1=0.1, frac2=0.5)
# 
# geweke.plot(model2.1$Sol, frac1=0.1,
#             frac2=0.5, nbins = 20,
#             pvalue = 0.05, auto.layout = TRUE)
# 
# 
# #geweke.diag(model2.1$VCV, frac1=0.1, frac2=0.5)
# 
# geweke.plot(model2.1$VCV, frac1=0.1,
#             frac2=0.5, nbins = 20,
#             pvalue = 0.05, auto.layout = TRUE)
# 
# 
# #raftery.diag(model2.1$Sol, q=0.05, r=0.01, s=0.95, converge.eps=0.001)
# 
# #raftery.diag(model2.1$VCV, q=0.05, r=0.01, s=0.95, converge.eps=0.001)
# 
# 
# # checking if there is autocorrelation between the values predicted
# 
# autocorr(model2.1$Sol)
# 
# autocorr(model2.1$VCV)
# 
# autocorr.plot(model2.1$Sol)
# 
# autocorr.plot(model2.1$VCV)
# 
# #heidel.diag(model2.1$Sol)
# #heidel.diag(model2.1$VCV)
# 
# # checking the results
# 
# # estimates for the fixed effects
# 
# posterior.mode(model2.1$Sol)
# 
# HPDinterval(model2.1$Sol)
# 
# 
# # variance explained by the pedigree and BirdID
# 
# posterior.mode(model2.1$VCV)
# 
# HPDinterval(model2.1$VCV)
# 
# 
# # estimating the percentage of variance explain by the pedigree and BirdID, i.e. checking Va and Vpe
# 
# model2.1.VP <- model2.1$VCV[,"animal"]+model2.1$VCV[,"BirdID"]+model2.1$VCV[,"units"]
# 
# model2.1.IDplusVA <- model2.1$VCV[,"animal"]+model2.1$VCV[,"BirdID"]
# 
# posterior.mode(model2.1.IDplusVA/model2.1.VP)
# 
# posterior.mode(model2.1$VCV[,"animal"]/model2.1.VP)
# 
# 
# # plotting the results as a forest plot
# 
# # list of values to be plotted
# 
# Modoff <- structure(c(posterior.mode(model2.1$VCV[,"animal"]/model2.1.VP), # mode heritability
#                       posterior.mode(model2.1.IDplusVA/model2.1.VP), # mode repeatability
#                       posterior.mode(model2.1$VCV[,"units"]/model2.1.VP), # mode residual
#                       HPDinterval(model2.1$VCV)[1,1]/posterior.mode(model2.1.VP), # lower CrI heritability
#                       (HPDinterval(model2.1$VCV)[1,1]+HPDinterval(model2.1$VCV)[2,1])/posterior.mode(model2.1.VP), # lower CrI repeatability
#                       HPDinterval(model2.1$VCV)[3,1]/posterior.mode(model2.1.VP), # lower CrI residuals
#                       HPDinterval(model2.1$VCV)[1,2]/posterior.mode(model2.1.VP), # upper CrI heritability
#                       (HPDinterval(model2.1$VCV)[1,2]+HPDinterval(model2.1$VCV)[2,2])/posterior.mode(model2.1.VP), # upper CrI repeatability
#                       HPDinterval(model2.1$VCV)[3,2]/posterior.mode(model2.1.VP)  # upper CrI residuals
# ), .Dim=c(3L, 3L), .Dimnames = list(c("heritability", "repeatability", "residuals"),
#                                     c("post.mode", "lower 95CrI", "upper 95CrI")))
# 
# 
# plot(1, 1,type="n",xlab="",ylab= "",xaxt="n",yaxt="n", bty="n")# using this nonsense because otherwise it overwrites existing plot
# 
# 
# 
# forestplot(c(" Va \n----- \n Vp","Va + Vpe \n-------------- \n      Vp"," Vr \n----- \n Vp"),
#            mean = Modoff[, "post.mode"], 
#            lower = Modoff[, "lower 95CrI"], 
#            upper = Modoff[, "upper 95CrI"],
#            boxsize = 0.05,
#            xticks=c(0, 0.1, 0.2, 0.3, 0.4,0.5,0.6,0.7,0.8,0.9,1.0), line.margin = unit(50,"mm"), 
#            lineheight =  unit(6, "cm"),
#            grid = TRUE,
#            lwd.xaxis = 2,lwd.zero = 2, lwd.ci = 2,
#            col=fpColors(box = "black", lines = "black",
#                         zero = "lightgray", text = "black", axes = "black",
#                         hrz_lines = "black"),
#            title="Standardized Elo-rating based on aggressive interactions\nfixed effects included: sex, age, tarsus\nVa=pedigree, Vpe=birdID, Vr=residuals, Vp=Va+Vpe+Vr")


######################################################################################################
# MODEL 3: StElo ~ SexEstimate + age + Tarsus.mean, random = ~ animal + BirdID + broodID
######################################################################################################

# # defining priors
# 
# prior3.1<-list(R = list(V = diag(1), nu = 0.002), 
#                G =list(G1 = list(V = diag(1), nu = 0.002), 
#                        G2 = list(V = diag(1), nu = 0.002),
#                        G3 = list(V = diag(1), nu = 0.002)))
# 
# 
# #  model
# 
# model3.1 <- MCMCglmm(StElo ~ 1,
#                      random = ~animal+BirdID+SocialBroodRef, 
#                      pedigree = prunedPed[,c(1,2,3)],
#                      data = h2dominance, 
#                      prior = prior3.1,
#                      nitt = 5000000,
#                      thin = 2500,
#                      burnin = 1250000,
#                      verbose=FALSE)
# 
# # save model if you haven't done it yet
# 
# save("model3.1", file = "models/mod1_h2_dominance_socialbrood.rda")


# loading the already run model

load("models/mod1_h2_dominance_socialbrood.rda")
summary(model3.1)


# evaluating how confident we can be that MCMCglmm found good answers

# checking how the model ran

plot(model3.1$Sol)
plot(model3.1$VCV)


# checking if the chain was more or less stable in the value

#geweke.diag(model3.1$Sol, frac1=0.1, frac2=0.5)

geweke.plot(model3.1$Sol, frac1=0.1,
            frac2=0.5, nbins = 20,
            pvalue = 0.05, auto.layout = TRUE)


#geweke.diag(model3.1$VCV, frac1=0.1, frac2=0.5)

geweke.plot(model3.1$VCV, frac1=0.1,
            frac2=0.5, nbins = 20,
            pvalue = 0.05, auto.layout = TRUE)


#raftery.diag(model3.1$Sol, q=0.05, r=0.01, s=0.95, converge.eps=0.001)

#raftery.diag(model3.1$VCV, q=0.05, r=0.01, s=0.95, converge.eps=0.001)


# checking if there is autocorrelation between the values predicted

autocorr(model3.1$Sol)

autocorr(model3.1$VCV)

autocorr.plot(model3.1$Sol)

autocorr.plot(model3.1$VCV)

#heidel.diag(model3.1$Sol)
#heidel.diag(model3.1$VCV)

# checking the results

# estimates for the fixed effects

posterior.mode(model3.1$Sol)

HPDinterval(model3.1$Sol)


# variance explained by the pedigree and BirdID

posterior.mode(model3.1$VCV)

HPDinterval(model3.1$VCV)


# estimating the percentage of variance explain by the pedigree and BirdID, i.e. checking Va and Vpe

model3.1.VP <- model3.1$VCV[,"animal"]+
  model3.1$VCV[,"BirdID"]+
  model3.1$VCV[,"SocialBroodRef"]+
  model3.1$VCV[,"units"]

# raw variance
posterior.mode(model3.1$VCV[,"animal"])
round(mean(model3.1$VCV[,"animal"]),0)
HPDinterval(model3.1$VCV[,"animal"])

posterior.mode(model3.1$VCV[,"BirdID"])
round(mean(model3.1$VCV[,"BirdID"]),0)
HPDinterval(model3.1$VCV[,"BirdID"])

posterior.mode(model3.1$VCV[,"SocialBroodRef"])
round(mean(model3.1$VCV[,"SocialBroodRef"]),0)
HPDinterval(model3.1$VCV[,"SocialBroodRef"])

posterior.mode(model3.1$VCV[,"units"])
round(mean(model3.1$VCV[,"units"]),0)
HPDinterval(model3.1$VCV[,"units"])

# proportion of variance
posterior.mode(model3.1$VCV[,"animal"]/model3.1.VP)
round(mean(model3.1$VCV[,"animal"]/model3.1.VP),3)
HPDinterval(model3.1$VCV[,"animal"]/model3.1.VP)

posterior.mode(model3.1$VCV[,"BirdID"]/model3.1.VP)
round(mean(model3.1$VCV[,"BirdID"]/model3.1.VP),3)
HPDinterval(model3.1$VCV[,"BirdID"]/model3.1.VP)

posterior.mode(model3.1$VCV[,"SocialBroodRef"]/model3.1.VP)
round(mean(model3.1$VCV[,"SocialBroodRef"]/model3.1.VP),3)
HPDinterval(model3.1$VCV[,"SocialBroodRef"]/model3.1.VP)

posterior.mode(model3.1$VCV[,"units"]/model3.1.VP)
round(mean(model3.1$VCV[,"units"]/model3.1.VP),3)
HPDinterval(model3.1$VCV[,"units"]/model3.1.VP)



# # plotting the results as a forest plot
# 
# # list of values to be plotted
# 
# Modoff <- structure(c(posterior.mode(model3.1$VCV[,"animal"]/model3.1.VP), # mode heritability
#                       posterior.mode(model3.1.IDplusVA/model3.1.VP), # mode repeatability
#                       posterior.mode(model3.1$VCV[,"SocialBroodRef"]/model3.1.VP), # mode brood
#                       posterior.mode(model3.1$VCV[,"units"]/model3.1.VP), # mode residual
#                       HPDinterval(model3.1$VCV)[1,1]/posterior.mode(model3.1.VP), # lower CrI heritability
#                       (HPDinterval(model3.1$VCV)[1,1]+HPDinterval(model3.1$VCV)[2,1])/posterior.mode(model3.1.VP), # lower CrI repeatability
#                       HPDinterval(model3.1$VCV)[3,1]/posterior.mode(model3.1.VP), # lower CrI residuals
#                       HPDinterval(model3.1$VCV)[4,1]/posterior.mode(model3.1.VP), # lower CrI brood
#                       HPDinterval(model3.1$VCV)[1,2]/posterior.mode(model3.1.VP), # upper CrI heritability
#                       (HPDinterval(model3.1$VCV)[1,2]+HPDinterval(model3.1$VCV)[2,2])/posterior.mode(model3.1.VP), # upper CrI repeatability
#                       HPDinterval(model3.1$VCV)[3,2]/posterior.mode(model3.1.VP), # upper CrI brood
#                       HPDinterval(model3.1$VCV)[4,2]/posterior.mode(model3.1.VP)# upper CrI residuals
#                       
# ), .Dim=c(4L, 3L), .Dimnames = list(c("heritability", "repeatability", "social brood" , "residuals"),
#                                     c("post.mode", "lower 95CrI", "upper 95CrI")))
# 
# 
# plot(1, 1,type="n",xlab="",ylab= "",xaxt="n",yaxt="n", bty="n")# using this nonsense because otherwise it overwrites existing plot
# 
# 
# 
# forestplot(c(" Va \n----- \n Vp","Va + Vpe \n-------------- \n      Vp"," Vs \n----- \n Vp"," Vr \n----- \n Vp"),
#            mean = Modoff[, "post.mode"], 
#            lower = Modoff[, "lower 95CrI"], 
#            upper = Modoff[, "upper 95CrI"],
#            boxsize = 0.05,
#            xticks=c(0, 0.1, 0.2, 0.3, 0.4,0.5,0.6,0.7,0.8,0.9,1.0), line.margin = unit(50,"mm"), 
#            lineheight =  unit(3.5, "cm"),
#            grid = TRUE,
#            lwd.xaxis = 2,lwd.zero = 2, lwd.ci = 2,
#            col=fpColors(box = "black", lines = "black",
#                         zero = "lightgray", text = "black", axes = "black",
#                         hrz_lines = "black"),
#            title="Standardized Elo-rating based on aggressive interactions\nfixed effects included: sex, age, tarsus\nVa=pedigree, Vpe=birdID, Vs=brood, Vr=residuals, Vp=Va+Vpe+Vs+Vr")



######################################################################################################
# MODEL 4: StElo ~ SexEstimate + age + Tarsus.mean, random = ~ animal + BirdID + socialMum
######################################################################################################

# defining priors

prior4.1<-list(R = list(V = diag(1), nu = 0.002), 
               G =list(G1 = list(V = diag(1), nu = 0.002), 
                       G2 = list(V = diag(1), nu = 0.002),
                       G3 = list(V = diag(1), nu = 0.002)))


# h2dominance.2 <- h2dominance[!(is.na(h2dominance$sex)) &
#                                !(is.na(h2dominance$age)) &
#                                !(is.na(h2dominance$tarsus)) &
#                                #!(is.na(h2dominance$dam)) &
#                                h2dominance$inped==1,]


# #  model
# 
# model4.1 <- MCMCglmm(StElo ~ 1,
#                      random = ~animal+BirdID+dam, 
#                      pedigree = prunedPed[,c(1,2,3)],
#                      data = h2dominance, 
#                      prior = prior3.1,
#                      nitt = 5000000,
#                      thin = 2500,
#                      burnin = 1250000,
#                      verbose=FALSE)
# 
# # save model if you haven't done it yet
# 
# save("model4.1", file = "models/mod1_h2_dominance_maternal.rda")


# loading the already run model

load("models/mod1_h2_dominance_maternal.rda")
summary(model4.1)


# evaluating how confident we can be that MCMCglmm found good answers

# checking how the model ran

plot(model4.1$Sol)
plot(model4.1$VCV)


# checking if the chain was more or less stable in the value

#geweke.diag(model4.1$Sol, frac1=0.1, frac2=0.5)

geweke.plot(model4.1$Sol, frac1=0.1,
            frac2=0.5, nbins = 20,
            pvalue = 0.05, auto.layout = TRUE)


#geweke.diag(model4.1$VCV, frac1=0.1, frac2=0.5)

geweke.plot(model4.1$VCV, frac1=0.1,
            frac2=0.5, nbins = 20,
            pvalue = 0.05, auto.layout = TRUE)


#raftery.diag(model4.1$Sol, q=0.05, r=0.01, s=0.95, converge.eps=0.001)

#raftery.diag(model4.1$VCV, q=0.05, r=0.01, s=0.95, converge.eps=0.001)


# checking if there is autocorrelation between the values predicted

autocorr(model4.1$Sol)

autocorr(model4.1$VCV)

autocorr.plot(model4.1$Sol)

autocorr.plot(model4.1$VCV)

#heidel.diag(model4.1$Sol)
#heidel.diag(model4.1$VCV)

# checking the results

# estimates for the fixed effects

posterior.mode(model4.1$Sol)

HPDinterval(model4.1$Sol)


# variance explained by the pedigree and BirdID

posterior.mode(model4.1$VCV)

HPDinterval(model4.1$VCV)


# estimating the percentage of variance explain by the pedigree and BirdID, i.e. checking Va and Vpe

model4.1.VP <- model4.1$VCV[,"animal"]+
  model4.1$VCV[,"BirdID"]+
  model4.1$VCV[,"dam"]+
  model4.1$VCV[,"units"]

# raw variance
posterior.mode(model4.1$VCV[,"animal"])
round(mean(model4.1$VCV[,"animal"]),0)
HPDinterval(model4.1$VCV[,"animal"])

posterior.mode(model4.1$VCV[,"BirdID"])
round(mean(model4.1$VCV[,"BirdID"]),0)
HPDinterval(model4.1$VCV[,"BirdID"])

posterior.mode(model4.1$VCV[,"dam"])
round(mean(model4.1$VCV[,"dam"]),0)
HPDinterval(model4.1$VCV[,"dam"])

posterior.mode(model4.1$VCV[,"units"])
round(mean(model4.1$VCV[,"units"]),0)
HPDinterval(model4.1$VCV[,"units"])

# proportion of variance
posterior.mode(model4.1$VCV[,"animal"]/model4.1.VP)
round(mean(model4.1$VCV[,"animal"]/model4.1.VP),3)
HPDinterval(model4.1$VCV[,"animal"]/model4.1.VP)

posterior.mode(model4.1$VCV[,"BirdID"]/model4.1.VP)
round(mean(model4.1$VCV[,"BirdID"]/model4.1.VP),3)
HPDinterval(model4.1$VCV[,"BirdID"]/model4.1.VP)

posterior.mode(model4.1$VCV[,"dam"]/model4.1.VP)
round(mean(model4.1$VCV[,"dam"]/model4.1.VP),3)
HPDinterval(model4.1$VCV[,"dam"]/model4.1.VP)

posterior.mode(model4.1$VCV[,"units"]/model4.1.VP)
round(mean(model4.1$VCV[,"units"]/model4.1.VP),3)
HPDinterval(model4.1$VCV[,"units"]/model4.1.VP)

# # plotting the results as a forest plot
# 
# # list of values to be plotted
# 
# Modoff <- structure(c(posterior.mode(model4.1$VCV[,"animal"]/model4.1.VP), # mode heritability
#                       posterior.mode(model4.1.IDplusVA/model4.1.VP), # mode repeatability
#                       posterior.mode(model4.1$VCV[,"dam"]/model4.1.VP), # mode brood
#                       posterior.mode(model4.1$VCV[,"units"]/model4.1.VP), # mode residual
#                       HPDinterval(model4.1$VCV)[1,1]/posterior.mode(model4.1.VP), # lower CrI heritability
#                       (HPDinterval(model4.1$VCV)[1,1]+HPDinterval(model4.1$VCV)[2,1])/posterior.mode(model4.1.VP), # lower CrI repeatability
#                       HPDinterval(model4.1$VCV)[3,1]/posterior.mode(model4.1.VP), # lower CrI residuals
#                       HPDinterval(model4.1$VCV)[4,1]/posterior.mode(model4.1.VP), # lower CrI brood
#                       HPDinterval(model4.1$VCV)[1,2]/posterior.mode(model4.1.VP), # upper CrI heritability
#                       (HPDinterval(model4.1$VCV)[1,2]+HPDinterval(model4.1$VCV)[2,2])/posterior.mode(model4.1.VP), # upper CrI repeatability
#                       HPDinterval(model4.1$VCV)[3,2]/posterior.mode(model4.1.VP), # upper CrI brood
#                       HPDinterval(model4.1$VCV)[4,2]/posterior.mode(model4.1.VP)# upper CrI residuals
#                       
# ), .Dim=c(4L, 3L), .Dimnames = list(c("heritability", "repeatability", "dam" , "residuals"),
#                                     c("post.mode", "lower 95CrI", "upper 95CrI")))
# 
# 
# plot(1, 1,type="n",xlab="",ylab= "",xaxt="n",yaxt="n", bty="n")# using this nonsense because otherwise it overwrites existing plot
# 
# 
# 
# forestplot(c(" Va \n----- \n Vp","Va + Vpe \n-------------- \n      Vp"," Vsm \n----- \n Vp"," Vr \n----- \n Vp"),
#            mean = Modoff[, "post.mode"], 
#            lower = Modoff[, "lower 95CrI"], 
#            upper = Modoff[, "upper 95CrI"],
#            boxsize = 0.05,
#            xticks=c(0, 0.1, 0.2, 0.3, 0.4,0.5,0.6,0.7,0.8,0.9,1.0), line.margin = unit(50,"mm"), 
#            lineheight =  unit(3.5, "cm"),
#            grid = TRUE,
#            lwd.xaxis = 2,lwd.zero = 2, lwd.ci = 2,
#            col=fpColors(box = "black", lines = "black",
#                         zero = "lightgray", text = "black", axes = "black",
#                         hrz_lines = "black"),
#            title="Standardized Elo-rating based on aggressive interactions\nfixed effects included: sex, age, tarsus\nVa=pedigree, Vpe=birdID, Vsm=socialMum, Vr=residuals, Vp=Va+Vpe+Vsm+Vr")
# 
# 

#########################################################
# plotting posterior distributions of heritability and PE
#########################################################

tiff("posterior_distributions_supplements.tiff", 
     height=15, width=30,
     units='cm', compression="lzw", res=600)

m <- rbind(c(1,2))


layout(m)

op <- par(oma = c(4,4,0.5,0.5) + 0.1,
          mar = c(1.5,3,1,0) + 0.1)

d.socialbrood <- density(model3.1$VCV[,"SocialBroodRef"]/model3.1.VP)

plot(d.socialbrood, 
     #type="n",
     xlab="",
     ylab= "",
     main="",
     xaxt="n",yaxt="n",xlim=c(0,0.2),ylim=c(0,800),
     family="serif")

axis(1,at=seq(0,0.2,0.05),
     las=1,
     cex.axis=1.5,
     family="serif") 

axis(2,at=seq(0,800,100),
     cex.axis=1.5,
     las=2,
     family="serif")

title(ylab = "density",
      line=1,
      outer = TRUE, cex.lab=3.25)

text(0.5,800,"(a)")

#polygon(d.socialbrood, col="grey75", border="black")


d.dam <- density(model4.1$VCV[,"dam"]/model4.1.VP)

plot(d.dam, 
     #type="n",
     main="",
     xlab="",
     ylab= "",
     xaxt="n",yaxt="n",xlim=c(0,0.2),ylim=c(0,800),
     family="serif")

axis(1,at=seq(0,0.2,0.05),
     las=1,
     cex.axis=1.5,
     family="serif") 

axis(2,at=seq(0,800,100),
     cex.axis=1.5,
     las=2,
     family="serif")

#polygon(d.dam, col="grey75", border="black")


title(xlab = "proportion of variance",
      outer = TRUE,
      line=2,
      cex.lab=3.25)

text(0.5,800,"(b)")

dev.off()


# ######################################################################################################
# # MODEL 5: StElo ~ SexEstimate + age + Tarsus.mean, random = ~ animal + BirdID + socialDad
# ######################################################################################################
# 
# # checking sample sizes
# 
# length(Data.v11$BirdID)
# length(unique(Data.v11$BirdID))
# 
# 
# # defining priors
# 
# prior5.1<-list(R = list(V = diag(1), nu = 0.002), 
#                G =list(G1 = list(V = diag(1), nu = 0.002), 
#                        G2 = list(V = diag(1), nu = 0.002),
#                        G3 = list(V = diag(1), nu = 0.002)))
# 
# 
# #  model
# 
# model5.1 <- MCMCglmm(StElo ~ scale(SexEstimate,scale=FALSE)+
#                        scale(age,scale=FALSE)+
#                        scale(Tarsus.mean,scale=FALSE), 
#                      random = ~animal+BirdID+SocialDadID, 
#                      pedigree = Ped,data = Data.v11, prior = prior5.1,
#                      #nitt = 1000, thin = 10, burnin = 100, verbose=FALSE)
#                      nitt = 4000000, thin = 3600, burnin = 400000, verbose=FALSE)
# 
# 
# # save model if you haven't done it yet
# 
# #save("model5.1", file = "mod5.1_StElo_h2andrepeatandsocialDad_fixedeffects_aggression_4million_3600thin-20160707.rda")
# 
# 
# # loading the already run model
# 
# load("mod5_mod5.1_StElo_h2andrepeatandsocialDad_fixedeffects_aggression_4million_3600thin-20160707.rda")
# summary(model5.1)
# 
# 
# # evaluating how confident we can be that MCMCglmm found good answers
# 
# # checking how the model ran
# 
# plot(model5.1$Sol)
# plot(model5.1$VCV)
# 
# 
# # checking if the chain was more or less stable in the value
# 
# #geweke.diag(model5.1$Sol, frac1=0.1, frac2=0.5)
# 
# geweke.plot(model5.1$Sol, frac1=0.1,
#             frac2=0.5, nbins = 20,
#             pvalue = 0.05, auto.layout = TRUE)
# 
# 
# #geweke.diag(model5.1$VCV, frac1=0.1, frac2=0.5)
# 
# geweke.plot(model5.1$VCV, frac1=0.1,
#             frac2=0.5, nbins = 20,
#             pvalue = 0.05, auto.layout = TRUE)
# 
# 
# #raftery.diag(model5.1$Sol, q=0.05, r=0.01, s=0.95, converge.eps=0.001)
# 
# #raftery.diag(model5.1$VCV, q=0.05, r=0.01, s=0.95, converge.eps=0.001)
# 
# 
# # checking if there is autocorrelation between the values predicted
# 
# autocorr(model5.1$Sol)
# 
# autocorr(model5.1$VCV)
# 
# autocorr.plot(model5.1$Sol)
# 
# autocorr.plot(model5.1$VCV)
# 
# #heidel.diag(model5.1$Sol)
# #heidel.diag(model5.1$VCV)
# 
# # checking the results
# 
# # estimates for the fixed effects
# 
# posterior.mode(model5.1$Sol)
# 
# HPDinterval(model5.1$Sol)
# 
# 
# # variance explained by the pedigree and BirdID
# 
# posterior.mode(model5.1$VCV)
# 
# HPDinterval(model5.1$VCV)
# 
# 
# # estimating the percentage of variance explain by the pedigree and BirdID, i.e. checking Va and Vpe
# 
# model5.1.VP <- model5.1$VCV[,"animal"]+model5.1$VCV[,"BirdID"]+model5.1$VCV[,"SocialDadID"]+model5.1$VCV[,"units"]
# 
# model5.1.IDplusVA <- model5.1$VCV[,"animal"]+model5.1$VCV[,"BirdID"]
# 
# posterior.mode(model5.1.IDplusVA/model5.1.VP)
# 
# posterior.mode(model5.1$VCV[,"animal"]/model5.1.VP)
# 
# posterior.mode(model5.1$VCV[,"SocialDadID"]/model5.1.VP)
# 
# 
# # plotting the results as a forest plot
# 
# # list of values to be plotted
# 
# Modoff <- structure(c(posterior.mode(model5.1$VCV[,"animal"]/model5.1.VP), # mode heritability
#                       posterior.mode(model5.1.IDplusVA/model5.1.VP), # mode repeatability
#                       posterior.mode(model5.1$VCV[,"SocialDadID"]/model5.1.VP), # mode brood
#                       posterior.mode(model5.1$VCV[,"units"]/model5.1.VP), # mode residual
#                       HPDinterval(model5.1$VCV)[1,1]/posterior.mode(model5.1.VP), # lower CrI heritability
#                       (HPDinterval(model5.1$VCV)[1,1]+HPDinterval(model5.1$VCV)[2,1])/posterior.mode(model5.1.VP), # lower CrI repeatability
#                       HPDinterval(model5.1$VCV)[3,1]/posterior.mode(model5.1.VP), # lower CrI residuals
#                       HPDinterval(model5.1$VCV)[4,1]/posterior.mode(model5.1.VP), # lower CrI brood
#                       HPDinterval(model5.1$VCV)[1,2]/posterior.mode(model5.1.VP), # upper CrI heritability
#                       (HPDinterval(model5.1$VCV)[1,2]+HPDinterval(model5.1$VCV)[2,2])/posterior.mode(model5.1.VP), # upper CrI repeatability
#                       HPDinterval(model5.1$VCV)[3,2]/posterior.mode(model5.1.VP), # upper CrI brood
#                       HPDinterval(model5.1$VCV)[4,2]/posterior.mode(model5.1.VP)# upper CrI residuals
#                       
# ), .Dim=c(4L, 3L), .Dimnames = list(c("heritability", "repeatability", "social Dad" , "residuals"),
#                                     c("post.mode", "lower 95CrI", "upper 95CrI")))
# 
# 
# plot(1, 1,type="n",xlab="",ylab= "",xaxt="n",yaxt="n", bty="n")# using this nonsense because otherwise it overwrites existing plot
# 
# 
# 
# forestplot(c(" Va \n----- \n Vp","Va + Vpe \n-------------- \n      Vp"," Vsd \n----- \n Vp"," Vr \n----- \n Vp"),
#            mean = Modoff[, "post.mode"], 
#            lower = Modoff[, "lower 95CrI"], 
#            upper = Modoff[, "upper 95CrI"],
#            boxsize = 0.05,
#            xticks=c(0, 0.1, 0.2, 0.3, 0.4,0.5,0.6,0.7,0.8,0.9,1.0), line.margin = unit(50,"mm"), 
#            lineheight =  unit(3.5, "cm"),
#            grid = TRUE,
#            lwd.xaxis = 2,lwd.zero = 2, lwd.ci = 2,
#            col=fpColors(box = "black", lines = "black",
#                         zero = "lightgray", text = "black", axes = "black",
#                         hrz_lines = "black"),
#            title="Standardized Elo-rating based on aggressive interactions\nfixed effects included: sex, age, tarsus\nVa=pedigree, Vpe=birdID, Vsd=socialDad, Vr=residuals, Vp=Va+Vpe+Vsd+Vr")
# 
# 
# 
# ######################################################################################################
# # MODEL 6: StElo ~ SexEstimate + age + Tarsus.mean, random = ~ animal + BirdID + GeneticMum
# ######################################################################################################
# 
# # checking sample sizes
# 
# length(Data.v12$BirdID)
# length(unique(Data.v12$BirdID))
# 
# 
# # defining priors
# 
# prior6.1<-list(R = list(V = diag(1), nu = 0.002), 
#                G =list(G1 = list(V = diag(1), nu = 0.002), 
#                        G2 = list(V = diag(1), nu = 0.002),
#                        G3 = list(V = diag(1), nu = 0.002)))
# 
# 
# #  model
# 
# model6.1 <- MCMCglmm(StElo ~ scale(SexEstimate,scale=FALSE)+
#                        scale(age,scale=FALSE)+
#                        scale(Tarsus.mean,scale=FALSE), 
#                      random = ~animal+BirdID+dam, 
#                      pedigree = Ped,data = Data.v12, prior = prior6.1,
#                      nitt = 1000, thin = 10, burnin = 100, verbose=FALSE)
# #nitt = 4000000, thin = 3600, burnin = 400000, verbose=FALSE)
# 
# 
# # save model if you haven't done it yet
# 
# #save("model6.1", file = "mod6.1_StElo_h2andrepeatandgeneticMum_fixedeffects_aggression_4million_3600thin-20160707.rda")
# 
# 
# # loading the already run model
# 
# load("mod6_mod6.1_StElo_h2andrepeatandgeneticMum_fixedeffects_aggression_4million_3600thin-20160707.rda")
# summary(model6.1)
# 
# 
# # evaluating how confident we can be that MCMCglmm found good answers
# 
# # checking how the model ran
# 
# plot(model6.1$Sol)
# plot(model6.1$VCV)
# 
# 
# # checking if the chain was more or less stable in the value
# 
# #geweke.diag(model6.1$Sol, frac1=0.1, frac2=0.5)
# 
# geweke.plot(model6.1$Sol, frac1=0.1,
#             frac2=0.5, nbins = 20,
#             pvalue = 0.05, auto.layout = TRUE)
# 
# 
# #geweke.diag(model6.1$VCV, frac1=0.1, frac2=0.5)
# 
# geweke.plot(model6.1$VCV, frac1=0.1,
#             frac2=0.5, nbins = 20,
#             pvalue = 0.05, auto.layout = TRUE)
# 
# 
# #raftery.diag(model6.1$Sol, q=0.05, r=0.01, s=0.95, converge.eps=0.001)
# 
# #raftery.diag(model6.1$VCV, q=0.05, r=0.01, s=0.95, converge.eps=0.001)
# 
# 
# # checking if there is autocorrelation between the values predicted
# 
# autocorr(model6.1$Sol)
# 
# autocorr(model6.1$VCV)
# 
# autocorr.plot(model6.1$Sol)
# 
# autocorr.plot(model6.1$VCV)
# 
# #heidel.diag(model6.1$Sol)
# #heidel.diag(model6.1$VCV)
# 
# # checking the results
# 
# # estimates for the fixed effects
# 
# posterior.mode(model6.1$Sol)
# 
# HPDinterval(model6.1$Sol)
# 
# 
# # variance explained by the pedigree and BirdID
# 
# posterior.mode(model6.1$VCV)
# 
# HPDinterval(model6.1$VCV)
# 
# 
# # estimating the percentage of variance explain by the pedigree and BirdID, i.e. checking Va and Vpe
# 
# model6.1.VP <- model6.1$VCV[,"animal"]+model6.1$VCV[,"BirdID"]+model6.1$VCV[,"dam"]+model6.1$VCV[,"units"]
# 
# model6.1.IDplusVA <- model6.1$VCV[,"animal"]+model6.1$VCV[,"BirdID"]
# 
# posterior.mode(model6.1.IDplusVA/model6.1.VP)
# 
# posterior.mode(model6.1$VCV[,"animal"]/model6.1.VP)
# 
# posterior.mode(model6.1$VCV[,"dam"]/model6.1.VP)
# 
# 
# # plotting the results as a forest plot
# 
# # list of values to be plotted
# 
# Modoff <- structure(c(posterior.mode(model6.1$VCV[,"animal"]/model6.1.VP), # mode heritability
#                       posterior.mode(model6.1.IDplusVA/model6.1.VP), # mode repeatability
#                       posterior.mode(model6.1$VCV[,"dam"]/model6.1.VP), # mode brood
#                       posterior.mode(model6.1$VCV[,"units"]/model6.1.VP), # mode residual
#                       HPDinterval(model6.1$VCV)[1,1]/posterior.mode(model6.1.VP), # lower CrI heritability
#                       (HPDinterval(model6.1$VCV)[1,1]+HPDinterval(model6.1$VCV)[2,1])/posterior.mode(model6.1.VP), # lower CrI repeatability
#                       HPDinterval(model6.1$VCV)[3,1]/posterior.mode(model6.1.VP), # lower CrI residuals
#                       HPDinterval(model6.1$VCV)[4,1]/posterior.mode(model6.1.VP), # lower CrI brood
#                       HPDinterval(model6.1$VCV)[1,2]/posterior.mode(model6.1.VP), # upper CrI heritability
#                       (HPDinterval(model6.1$VCV)[1,2]+HPDinterval(model6.1$VCV)[2,2])/posterior.mode(model6.1.VP), # upper CrI repeatability
#                       HPDinterval(model6.1$VCV)[3,2]/posterior.mode(model6.1.VP), # upper CrI brood
#                       HPDinterval(model6.1$VCV)[4,2]/posterior.mode(model6.1.VP)# upper CrI residuals
#                       
# ), .Dim=c(4L, 3L), .Dimnames = list(c("heritability", "repeatability", "genetic Mum" , "residuals"),
#                                     c("post.mode", "lower 95CrI", "upper 95CrI")))
# 
# 
# plot(1, 1,type="n",xlab="",ylab= "",xaxt="n",yaxt="n", bty="n")# using this nonsense because otherwise it overwrites existing plot
# 
# 
# 
# forestplot(c(" Va \n----- \n Vp","Va + Vpe \n-------------- \n      Vp"," Vgm \n----- \n Vp"," Vr \n----- \n Vp"),
#            mean = Modoff[, "post.mode"], 
#            lower = Modoff[, "lower 95CrI"], 
#            upper = Modoff[, "upper 95CrI"],
#            boxsize = 0.05,
#            xticks=c(0, 0.1, 0.2, 0.3, 0.4,0.5,0.6,0.7,0.8,0.9,1.0), line.margin = unit(50,"mm"), 
#            lineheight =  unit(3.5, "cm"),
#            grid = TRUE,
#            lwd.xaxis = 2,lwd.zero = 2, lwd.ci = 2,
#            col=fpColors(box = "black", lines = "black",
#                         zero = "lightgray", text = "black", axes = "black",
#                         hrz_lines = "black"),
#            title="Standardized Elo-rating based on aggressive interactions\nfixed effects included: sex, age, tarsus\nVa=pedigree, Vpe=birdID, Vgm=geneticMum, Vr=residuals, Vp=Va+Vpe+Vgm+Vr")
# 
# 
# ######################################################################################################
# # MODEL 7: StElo ~ SexEstimate + age + Tarsus.mean, random = ~ animal + BirdID + Nestbox
# ######################################################################################################
# 
# # checking sample sizes
# 
# length(Data.v13$BirdID)
# length(unique(Data.v13$BirdID))
# 
# 
# # defining priors
# 
# prior7.1<-list(R = list(V = diag(1), nu = 0.002), 
#                G =list(G1 = list(V = diag(1), nu = 0.002), 
#                        G2 = list(V = diag(1), nu = 0.002),
#                        G3 = list(V = diag(1), nu = 0.002)))
# 
# 
# #  model
# 
# model7.1 <- MCMCglmm(StElo ~ scale(SexEstimate,scale=FALSE)+
#                        scale(age,scale=FALSE)+
#                        scale(Tarsus.mean,scale=FALSE), 
#                      random = ~animal+BirdID+NestboxName, 
#                      pedigree = Ped,data = Data.v13, prior = prior7.1,
#                      #nitt = 1000, thin = 10, burnin = 100, verbose=FALSE)
#                      nitt = 4000000, thin = 3600, burnin = 400000, verbose=FALSE)
# 
# 
# # save model if you haven't done it yet
# 
# save("model7.1", file = "mod7_mod7.1_StElo_h2andrepeatandNestbox_fixedeffects_aggression_4million_3600thin-20160708.rda")
# 
# 
# # loading the already run model
# 
# #load("mod7_mod7.1_StElo_h2andrepeatandNestbox_fixedeffects_aggression_4million_3600thin-20160708.rda")
# #summary(model7.1)
# 
# 
# # evaluating how confident we can be that MCMCglmm found good answers
# 
# # checking how the model ran
# 
# plot(model7.1$Sol)
# plot(model7.1$VCV)
# 
# 
# # checking if the chain was more or less stable in the value
# 
# #geweke.diag(model7.1$Sol, frac1=0.1, frac2=0.5)
# 
# geweke.plot(model7.1$Sol, frac1=0.1,
#             frac2=0.5, nbins = 20,
#             pvalue = 0.05, auto.layout = TRUE)
# 
# 
# #geweke.diag(model7.1$VCV, frac1=0.1, frac2=0.5)
# 
# geweke.plot(model7.1$VCV, frac1=0.1,
#             frac2=0.5, nbins = 20,
#             pvalue = 0.05, auto.layout = TRUE)
# 
# 
# #raftery.diag(model7.1$Sol, q=0.05, r=0.01, s=0.95, converge.eps=0.001)
# 
# #raftery.diag(model7.1$VCV, q=0.05, r=0.01, s=0.95, converge.eps=0.001)
# 
# 
# # checking if there is autocorrelation between the values predicted
# 
# autocorr(model7.1$Sol)
# 
# autocorr(model7.1$VCV)
# 
# autocorr.plot(model7.1$Sol)
# 
# autocorr.plot(model7.1$VCV)
# 
# #heidel.diag(model7.1$Sol)
# #heidel.diag(model7.1$VCV)
# 
# # checking the results
# 
# # estimates for the fixed effects
# 
# posterior.mode(model7.1$Sol)
# 
# HPDinterval(model7.1$Sol)
# 
# 
# # variance explained by the pedigree and BirdID
# 
# posterior.mode(model7.1$VCV)
# 
# HPDinterval(model7.1$VCV)
# 
# 
# # estimating the percentage of variance explain by the pedigree and BirdID, i.e. checking Va and Vpe
# 
# model7.1.VP <- model7.1$VCV[,"animal"]+model7.1$VCV[,"BirdID"]+model7.1$VCV[,"NestboxName"]+model7.1$VCV[,"units"]
# 
# model7.1.IDplusVA <- model7.1$VCV[,"animal"]+model7.1$VCV[,"BirdID"]
# 
# posterior.mode(model7.1.IDplusVA/model7.1.VP)
# 
# posterior.mode(model7.1$VCV[,"animal"]/model7.1.VP)
# 
# posterior.mode(model7.1$VCV[,"NestboxName"]/model7.1.VP)
# 
# 
# # plotting the results as a forest plot
# 
# # list of values to be plotted
# 
# Modoff <- structure(c(posterior.mode(model7.1$VCV[,"animal"]/model7.1.VP), # mode heritability
#                       posterior.mode(model7.1.IDplusVA/model7.1.VP), # mode repeatability
#                       posterior.mode(model7.1$VCV[,"NestboxName"]/model7.1.VP), # mode NestboxName
#                       posterior.mode(model7.1$VCV[,"units"]/model7.1.VP), # mode residual
#                       HPDinterval(model7.1$VCV)[1,1]/posterior.mode(model7.1.VP), # lower CrI heritability
#                       (HPDinterval(model7.1$VCV)[1,1]+HPDinterval(model7.1$VCV)[2,1])/posterior.mode(model7.1.VP), # lower CrI repeatability
#                       HPDinterval(model7.1$VCV)[3,1]/posterior.mode(model7.1.VP), # lower CrI residuals
#                       HPDinterval(model7.1$VCV)[4,1]/posterior.mode(model7.1.VP), # lower CrI NestboxName
#                       HPDinterval(model7.1$VCV)[1,2]/posterior.mode(model7.1.VP), # upper CrI heritability
#                       (HPDinterval(model7.1$VCV)[1,2]+HPDinterval(model7.1$VCV)[2,2])/posterior.mode(model7.1.VP), # upper CrI repeatability
#                       HPDinterval(model7.1$VCV)[3,2]/posterior.mode(model7.1.VP), # upper CrI NestboxName
#                       HPDinterval(model7.1$VCV)[4,2]/posterior.mode(model7.1.VP)# upper CrI residuals
#                       
# ), .Dim=c(4L, 3L), .Dimnames = list(c("heritability", "repeatability", "nestbox" , "residuals"),
#                                     c("post.mode", "lower 95CrI", "upper 95CrI")))
# 
# 
# plot(1, 1,type="n",xlab="",ylab= "",xaxt="n",yaxt="n", bty="n")# using this nonsense because otherwise it overwrites existing plot
# 
# 
# 
# forestplot(c(" Va \n----- \n Vp","Va + Vpe \n-------------- \n      Vp"," Vnb \n----- \n Vp"," Vr \n----- \n Vp"),
#            mean = Modoff[, "post.mode"], 
#            lower = Modoff[, "lower 95CrI"], 
#            upper = Modoff[, "upper 95CrI"],
#            boxsize = 0.05,
#            xticks=c(0, 0.1, 0.2, 0.3, 0.4,0.5,0.6,0.7,0.8,0.9,1.0), line.margin = unit(50,"mm"), 
#            lineheight =  unit(3.5, "cm"),
#            grid = TRUE,
#            lwd.xaxis = 2,lwd.zero = 2, lwd.ci = 2,
#            col=fpColors(box = "black", lines = "black",
#                         zero = "lightgray", text = "black", axes = "black",
#                         hrz_lines = "black"),
#            title="Standardized Elo-rating based on aggressive interactions\nfixed effects included: sex, age, tarsus\nVa=pedigree, Vpe=birdID, Vnb=nestbox, Vr=residuals, Vp=Va+Vpe+Vgm+Vr")
# 
# 
# ######################################################################################################
# # MODEL 8: StElo ~ SexEstimate + age + Tarsus.mean + season, random = ~ animal + BirdID
# ######################################################################################################
# 
# # checking sample sizes
# length(Data2.1$BirdID)
# length(unique(Data2.1$BirdID))
# 
# 
# # defining priors
# 
# prior8.1<-list(R = list(V = diag(1), nu = 0.002), 
#                G =list(G1 = list(V = diag(1), nu = 0.002), G2 = list(V = diag(1), nu = 0.002)))
# 
# 
# #  model
# 
# model8.1 <- MCMCglmm(StElo ~ scale(SexEstimate,scale=FALSE)+scale(age,scale=FALSE)+scale(Tarsus.mean,scale=FALSE)
#                      +scale(season,scale=FALSE), 
#                      random = ~animal+BirdID, pedigree = Ped,data = Data2.1, prior = prior8.1,
#                      nitt = 7500000, thin = 3600, burnin = 750000, verbose=FALSE)
# #nitt = 1000, thin = 10, burnin = 100, verbose=FALSE)
# 
# 
# # save model if you haven't done it yet
# 
# save("model8.1", file = "mod8_mod8.1_StElo_h2andrepeatand_fixedeffectsv2_aggression_7.5million_3600thin-20160709.rda")
# 
# 
# # loading the already run model
# 
# #load("mod8_mod8.1_StElo_h2andrepeatand_fixedeffectsv2_aggression_4million_3600thin-20160709.rda")
# #summary(model8.1)
# 
# 
# # evaluating how confident we can be that MCMCglmm found good answers
# 
# # checking how the model ran
# 
# plot(model8.1$Sol)
# plot(model8.1$VCV)
# 
# 
# # checking if the chain was more or less stable in the value
# 
# #geweke.diag(model8.1$Sol, frac1=0.1, frac2=0.5)
# 
# geweke.plot(model8.1$Sol, frac1=0.1,
#             frac2=0.5, nbins = 20,
#             pvalue = 0.05, auto.layout = TRUE)
# 
# 
# #geweke.diag(model8.1$VCV, frac1=0.1, frac2=0.5)
# 
# geweke.plot(model8.1$VCV, frac1=0.1,
#             frac2=0.5, nbins = 20,
#             pvalue = 0.05, auto.layout = TRUE)
# 
# 
# #raftery.diag(model8.1$Sol, q=0.05, r=0.01, s=0.95, converge.eps=0.001)
# 
# #raftery.diag(model8.1$VCV, q=0.05, r=0.01, s=0.95, converge.eps=0.001)
# 
# 
# # checking if there is autocorrelation between the values predicted
# 
# autocorr(model8.1$Sol)
# 
# autocorr(model8.1$VCV)
# 
# autocorr.plot(model8.1$Sol)
# 
# autocorr.plot(model8.1$VCV)
# 
# #heidel.diag(model8.1$Sol)
# #heidel.diag(model8.1$VCV)
# 
# # checking the results
# 
# # estimates for the fixed effects
# 
# posterior.mode(model8.1$Sol)
# 
# HPDinterval(model8.1$Sol)
# 
# 
# # variance explained by the pedigree and BirdID
# 
# posterior.mode(model8.1$VCV)
# 
# HPDinterval(model8.1$VCV)
# 
# 
# # estimating the percentage of variance explain by the pedigree and BirdID, i.e. checking Va and Vpe
# 
# model8.1.VP <- model8.1$VCV[,"animal"]+model8.1$VCV[,"BirdID"]+model8.1$VCV[,"units"]
# 
# model8.1.IDplusVA <- model8.1$VCV[,"animal"]+model8.1$VCV[,"BirdID"]
# 
# posterior.mode(model8.1.IDplusVA/model8.1.VP)
# 
# posterior.mode(model8.1$VCV[,"animal"]/model8.1.VP)
# 
# 
# # plotting the results as a forest plot
# 
# # list of values to be plotted
# 
# Modoff <- structure(c(posterior.mode(model8.1$VCV[,"animal"]/model8.1.VP), # mode heritability
#                       posterior.mode(model8.1.IDplusVA/model8.1.VP), # mode repeatability
#                       posterior.mode(model8.1$VCV[,"units"]/model8.1.VP), # mode residual
#                       HPDinterval(model8.1$VCV)[1,1]/posterior.mode(model8.1.VP), # lower CrI heritability
#                       (HPDinterval(model8.1$VCV)[1,1]+HPDinterval(model8.1$VCV)[2,1])/posterior.mode(model8.1.VP), # lower CrI repeatability
#                       HPDinterval(model8.1$VCV)[3,1]/posterior.mode(model8.1.VP), # lower CrI residuals
#                       HPDinterval(model8.1$VCV)[1,2]/posterior.mode(model8.1.VP), # upper CrI heritability
#                       (HPDinterval(model8.1$VCV)[1,2]+HPDinterval(model8.1$VCV)[2,2])/posterior.mode(model8.1.VP), # upper CrI repeatability
#                       HPDinterval(model8.1$VCV)[3,2]/posterior.mode(model8.1.VP) # upper CrI NestboxName
#                       
# ), .Dim=c(3L, 3L), .Dimnames = list(c("heritability", "repeatability", "residuals"),
#                                     c("post.mode", "lower 95CrI", "upper 95CrI")))
# 
# 
# plot(1, 1,type="n",xlab="",ylab= "",xaxt="n",yaxt="n", bty="n")# using this nonsense because otherwise it overwrites existing plot
# 
# 
# 
# forestplot(c(" Va \n----- \n Vp","Va + Vpe \n-------------- \n      Vp", "Vr \n----- \n Vp"),
#            mean = Modoff[, "post.mode"], 
#            lower = Modoff[, "lower 95CrI"], 
#            upper = Modoff[, "upper 95CrI"],
#            boxsize = 0.05,
#            xticks=c(0, 0.1, 0.2, 0.3, 0.4,0.5,0.6,0.7,0.8,0.9,1.0), line.margin = unit(50,"mm"), 
#            lineheight =  unit(3.5, "cm"),
#            grid = TRUE,
#            lwd.xaxis = 2,lwd.zero = 2, lwd.ci = 2,
#            col=fpColors(box = "black", lines = "black",
#                         zero = "lightgray", text = "black", axes = "black",
#                         hrz_lines = "black"),
#            title="Standardized Elo-rating based on aggressive interactions\nfixed effects included: sex, age, tarsus, season\nVa=pedigree, Vpe=birdID, Vr=residuals, Vp=Va+Vpe+Vr")
