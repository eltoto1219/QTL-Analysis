# R script for analyzing werg QTL data
# Antonio Mendoza
# werg.r script
# March 18, 2017
#
# Clean up and set directory
rm(list=ls())
getwd()
setwd("/Users/Antonio/Desktop/Gotwals/Intro/QTL")
#
#Making Functions
histCurve <- function(h) {
  myhist<- hist(as.matrix(h), probability = T, xlim = c(-4, 4), main = paste("Hisotgram of ",
                                                                             deparse(substitute(h)), sep =""))
  mycurve <- curve(dnorm(x, mean=mean(as.matrix(h), na.rm = TRUE),
                         sd=sd(as.matrix(h), na.rm = TRUE)),
                   y =-4,to=4, add=T)}
#
# Load the library and data
library("qtl")
werg <- read.cross("csv",
                   file="Werg.csv",
                   genotypes = c("N", "R", "H"),
                   na.strings = "-",
                   alleles=c("N", "R"))
#
#checking phenotypes
names(werg$pheno)
#
# Jittermap will move the genetic markers apart slightly so the results are better.
jittermap(werg)
#
# Look at the data, make sure it’s pretty clean (No red spots, especially in bottom right corner)
werg <- est.rf(werg)
jpeg( "jittermap.jpg")
plot.rf(werg)
graphics.off()
#
#Genetic map (the horizontal lines are genetic markers) & look at missing data
jpeg( "genetic_maps.jpg", width = 800, height = 400)
par(mfrow=c(1,2), pin = c(4,4))
plot.map(werg)
plot.missing(werg)
graphics.off()
dev.off()
#
#setting up histograms
Ftot<- werg$pheno$Ftotden
P<- werg$pheno$P
Crt<- werg$pheno$CrtThk
#Root mean squared scaling and mean-centering
Ftot_normalized <- scale(Ftot)
P_normalized <- scale(P)
Crt_normalized <- scale(Crt)
#hists
jpeg( "hists_raw_normalized.jpg", height = 480, width = 580)
par(mfrow = c(2, 3), pin = c(2,2))
hist(Ftot)
hist(P)
hist(Crt)
histCurve(Ftot_normalized)
histCurve(P_normalized)
histCurve(Crt_normalized)
graphics.off()
dev.off()
#
#QQ Plots
jpeg( "qqplots.jpg", width = 1200, height = 400)
par(mfrow = c(1, 3), pin = c(4,4))
qqnorm(Ftot_normalized, main = " QQ Plot of Ftot")
qqline(Ftot_normalized )
qqnorm(P_normalized, main = " QQ Plot ofP")
qqline(P_normalized)
qqnorm(Crt_normalized, main = " QQ Plot of Crt")
qqline(Crt_normalized)
qqline(Crt_normalized)
graphics.off()
dev.off()
#
#setting up for a mainscans
#Calculate genetic probability map to see what scan SHOULD look like, then run a simulation
werg <- calc.genoprob(werg, step=2.0, off.end=0.0, error.prob=1.0e-4, map.function="haldane",
                      stepwidth="fixed")
werg <- sim.geno(werg, n.draws=32, step=2.0, off.end=0.0, error.prob=1.0e-4,
                 map.function="haldane",
                 stepwidth="fixed")
#
# Mainscans for QTL & run for 100 permutations per each phenotype
werg_Ftot <- scanone(werg, pheno.col=3, model="normal", method="em")
werg_Ftot_perm <- scanone(werg, pheno.col=3, model="normal", method="em", n.perm=100)
werg_P <- scanone(werg, pheno.col=4, model="normal", method="em")
werg_P_perm <- scanone(werg, pheno.col=4, model="normal", method="em", n.perm=100)
werg_Crt <- scanone(werg, pheno.col=5, model="normal", method="em")
werg_Crt_perm <- scanone(werg, pheno.col=5, model="normal", method="em", n.perm=100)
#
# Plot the mainscans with threshold like 63%, 13%, and 5% confidences
jpeg( "mainscans.jpg", width = 1200, height = 400)
par(mfrow = c(1, 3), pin = c(4,4))
thresh <- summary(werg_Ftot_perm, alpha=c(0.63, 0.10, 0.05))
plot(werg_Ftot, main="Mainscan Plot of Ftot")
abline(h=thresh[1], col="blue")
abline(h=thresh[2], col="green")
abline(h=thresh[3], col="red")
#
thresh <- summary(werg_P_perm, alpha=c(0.63, 0.10, 0.05))
plot(werg_P, main="Mainscan Plot of P")
abline(h=thresh[1], col="blue")
abline(h=thresh[2], col="green")
abline(h=thresh[3], col="red")
#
thresh <- summary(werg_Crt_perm, alpha=c(0.63, 0.10, 0.05))
plot(werg_Crt, main="Mainscan Plot of Crt")
abline(h=thresh[1], col="blue")
abline(h=thresh[2], col="green")
abline(h=thresh[3], col="red")
graphics.off()
dev.off()
#
#Summary output of mainscan in text
summary(werg_Ftot, perm=werg_Ftot_perm, lodcolumn=1, alpha=0.05 )
summary(werg_P, perm=werg_P_perm, lodcolumn=1, alpha=0.05 )
summary(werg_Crt, perm=werg_Crt_perm, lodcolumn=1, alpha=0.05 )
#
#Setting up effect plots for higest LOD scores
jpeg( "effects.jpg", width = 1200, height = 400)
par(mfrow = c(1, 3), pin = c(4,4))
firsttdeffect <- find.marker(werg, chr=11, pos=44.7)
effectplot(werg, pheno.col=3, mname1=firsttdeffect, main = "Effect Plot Ftot")
seconddeffect <- find.marker(werg, chr=12, pos=5.52 )
effectplot(werg, pheno.col=4, mname1=seconddeffect, main = "Effect Plot of P")
thirdeffect <- find.marker(werg, chr=4, pos=69.05)
effectplot(werg, pheno.col=5, mname1=thirdeffect, main = "Effect Plot of Crt")
graphics.off()
Printed by Wolfram Mathematica Student Edition
graphics.off()
dev.off()
#