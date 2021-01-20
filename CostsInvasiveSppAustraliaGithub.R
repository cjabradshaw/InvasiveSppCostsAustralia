################################################################
## R Code to calculate costs of invasive species in Australia ##
## Corey Bradshaw, Andrew Hoskins, Phillip Haubrock, Boris Leroy
## corey.bradshaw@flinders.edu.au
## Flinders University, Australia
## January 2021
################################################################

## Accompanies paper: Bradshaw, CJA, AJ Hoskins, PJ Haubrock, RN Cuthbert, C Diagne, B Leroy, L Andrews, B Page,
##                      P Cassey, AW Sheppard, F Courchamp. In review. Detailed assessment of the reported economic costs of
##                      invasive alien species in Australia. NeoBiota


rm(list=ls())

# Set functions
AICc <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]
  }
  return(AICc.vec)
}

delta.AIC <- function(x) x - min(x) ## where x is a vector of AIC
weight.AIC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dAIC
ch.dev <- function(x) ((( as.numeric(x$null.deviance) - as.numeric(x$deviance) )/ as.numeric(x$null.deviance))*100) ## % change in deviance, where x is glm object

linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
  fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
  AIC.vec <- c(AICc(fit.full),AICc(fit.null))
  dAIC.vec <- delta.AIC(AIC.vec); wAIC.vec <- weight.AIC(dAIC.vec)
  ER <- wAIC.vec[1]/wAIC.vec[2]
  r.sq.adj <- as.numeric(summary(fit.full)[9])
  return(c(ER,r.sq.adj))
}

# set directory
library(sandwich)
library(ggalluvial)
library(lmtest)
library(robustbase)
library(ggplot2)
library(ggpubr)
library(plyr)
library(dplyr)
library(grid)
library(invacost) # obtained from Boris Leroy (leroy.boris@gmail.com)

## load pre-filtered file from InvaCost 2.1 & herbivore database combined (not available until paper is published)
invacost <- read.csv("AustraliaUpdate3.csv",sep=",", stringsAsFactors=TRUE)
dim(invacost)

# fix environment entries
invacost$Environment2 <- invacost$Environment 
invacost <- within(invacost, Environment2[Species == 'Alternanthera philoxeroides' & Environment2 == 'Aquatic/Terrestrial'] <- 'Aquatic') 
invacost <- within(invacost, Environment2[Species == 'Alternanthera philoxeroides' & Environment2 == 'Diverse/Unspecified'] <- 'Aquatic') 
invacost <- within(invacost, Environment2[Environment2 == 'Aquatic/Terrestrial'] <- 'Diverse/Unspecified') 
invacost <- within(invacost, Environment2[Species == 'Rhinella marina' & Environment2 == 'Terrestrial'] <- 'Aquatic') 
invacost <- within(invacost, Environment2[Common_name == 'Mammals/Amphibians/Insects/Weeds' & Environment2 == 'Terrestrial'] <- 'Diverse/Unspecified') 
invacost <- within(invacost, Environment2[Common_name == 'Pest animals vertebrates (including aquatic species)' & Environment2 == 'Terrestrial'] <- 'Diverse/Unspecified')

# remove avoided costs
invacost <- invacost[which(invacost$Type_of_cost != "Avoided costs"),]
dim(invacost)

invacost$Type_of_cost <- revalue(invacost$Type_of_cost, c("Control-cost" = "Management",
                                                "Control-costs" = "Management",
                                                "Contol/Management" = "Management",
                                                "Control/Mamagement" = "Management",
                                                "Contol/research" = "Management",
                                                "Control cost" = "Management",
                                                "Control" = "Management",
                                                "Management/Repair" = "Management",
                                                "Medical care" = "Management",
                                                "Export facilitation" = "Management",
                                                "direct control cost" = "Management",
                                                "Control/research" = "Management",
                                                "Control/Damage-Loss" = "Mixed",
                                                "Control/Education" = "Management",
                                                "Control/Medical care" = "Mixed",
                                                "Control/Prevention" = "Management",
                                                "Control/Monitoring" = "Management",
                                                "Control/Research" = "Management",
                                                "Control/Research/Surveillance" = "Management",
                                                "Control/Surveillance" = "Management",
                                                "Control/Management" = "Management",
                                                "Damage-Loss" = "Damage",
                                                "Damage-Loss/Management" = "Mixed",
                                                "Damage-loss, Control-costs" = "Mixed",
                                                "Damage-costs/Control-costs" = "Mixed",
                                                "Damage-costs" = "Damage",
                                                "Damage-loss" = "Damage",
                                                "Damage repair/Eradication" = "Mixed",
                                                "Diverse/Unspecified" = "Mixed",
                                                "Early detection/Management" = "Management",
                                                "Education" = "Management",
                                                "Eradication" = "Management",
                                                "Funding" = "Management",
                                                "Monitoring " = "Management",
                                                "Monitoring/Research" = "Management",
                                                "Personal loss" = "Damage",
                                                "Prevention" = "Management",
                                                "Prevention/Research" = "Management",
                                                "Prevention/Surveillance" = "Management",
                                                "Research" = "Management"))

# check
table(invacost$Type_of_cost)

invacost$Impacted_sector <- revalue(invacost$Impacted_sector, c("Agriculture/Authorities-Stakeholders" = "Mixed",
                                                      "Agriculture/Authorities-Stakeholders/Public and social welfare" = "Mixed",
                                                      "Agriculture/Environment" = "Mixed",
                                                      "Agricultural/forestry" = "Mixed",
                                                      "Agricultural" = "Agriculture",
                                                      "Agriculture/Health/Environment/Trade" = "Mixed",
                                                      "Agriculture/Trade" = "Mixed",
                                                      "Agriculture/Environment/Health/Public and social welfare" = "Mixed",
                                                      "Agriculture/Health/Environment" = "Mixed",
                                                      "Agriculture/Public and social welfare" = "Mixed",
                                                      "Authorities-Stakeholders/Health/Public and social welfare" = "Mixed",
                                                      "Authorities-Stakeholders/Public and social welfare" = "Mixed",
                                                      "cost of deer grazing pressure" = "Agriculture",
                                                      "culling operations by landowners" = "Agriculture",
                                                      "Diverse/Unspecified" = "Mixed",
                                                      "Public health" = "Health",
                                                      "Environment/Public and social welfare" = "Mixed",
                                                      "landholder survey" = "Public and social welfare",
                                                      "state estimate based on landholder survey" = "Public and social welfare",
                                                      "Unspecified" = "Mixed"))


# check
table(invacost$Impacted_sector)


#expanding the total costs for the invasives 
Australia1 <- expandYearlyCosts(invacost, startcolumn = "Probable_starting_year_low_margin", endcolumn = "Probable_ending_year_low_margin")
options(scipen=999)

# delete costs before 1960 and after 2017 
Australia2 <- Australia1 %>% filter(Impact_year <= "2020")
Australia <- Australia2 %>% filter(Impact_year >= "1960")

# Type_of_Cost
Australia$Type_of_cost <- revalue(Australia$Type_of_cost, c("Damage_costs" = "Damage-Loss",
                                                      "Mixed_costs"="Diverse/Unspecified"))
# Impacted_sector
Australia$Impacted_sector <- revalue(Australia$Impacted_sector, c("Authorities-Stakeholders/Public and social welfare" = "Mixed",
   
                                                                  "Unspecified"="Mixed"))
## cost calcs in US$ billion
Australia$cost <- as.numeric(gsub(",", "", Australia$Cost_estimate_per_year_2017_USD_exchange_rate))
Australia <- Australia[!is.na(Australia$cost),]
Australia$cost_bil <- (Australia$cost/1000000000)
sum(Australia$cost_bil) # in USD
1.304796*sum(Australia$cost_bil) # in 2017 AUD
length(Australia$cost_bil)

# implementation
AustraliaObs <- Australia[Australia$Implementation %in% c("Observed"),] 
sum(AustraliaObs$cost_bil)
length(AustraliaObs$cost_bil)
sum(AustraliaObs$cost_bil) / sum(Australia$cost_bil) # proportion

# filtering to potential costs only
AustraliaPot <- Australia[Australia$Implementation %in% c("Potential"),]
sum(AustraliaPot$cost_bil)
length(AustraliaPot$cost_bil)
sum(AustraliaPot$cost_bil) / sum(Australia$cost_bil) # proportion

#filtering to highly reliable (observed)
AustraliaHigh <- AustraliaObs[AustraliaObs$Method_reliability %in% c("High"),]
sum(AustraliaHigh$cost_bil)
sum(AustraliaHigh$cost_bil) / sum(Australia$cost_bil) # proportion
length(AustraliaHigh$cost_bil)

#filtering to low-reliability observed costs only
AustraliaLow <- AustraliaObs[AustraliaObs$Method_reliability %in% c("Low"),]
sum(AustraliaLow$cost_bil)
length(AustraliaLow$cost_bil)


## SUMMARIES
# environment
aggregate(cost_bil~Environment,data=Australia,FUN="sum")
aggregate(cost_bil~Environment,data=Australia,FUN="length")

100 * (aggregate(cost_bil~Environment,data=Australia,FUN="sum")$cost_bil / (sum(aggregate(cost_bil~Environment,data=Australia,FUN="sum")$cost_bil)) )

aggregate(cost_bil~Environment,data=AustraliaHigh,FUN="sum")
aggregate(cost_bil~Environment,data=AustraliaHigh,FUN="sum")$cost_bil / sum(aggregate(cost_bil~Environment,data=AustraliaHigh,FUN="sum")$cost_bil)
aggregate(cost_bil~Environment,data=AustraliaHigh,FUN="length")

# type
aggregate(cost_bil~Type_of_cost,data=Australia,FUN=sum)
sum(aggregate(cost_bil~Type_of_cost,data=Australia,FUN=sum)$cost_bil)
aggregate(cost_bil~Type_of_cost,data=Australia,FUN=length)

aggregate(cost_bil~Type_of_cost,data=AustraliaHigh,FUN=sum)
aggregate(cost_bil~Type_of_cost,data=AustraliaHigh,FUN=sum)$cost_bil / sum(aggregate(cost_bil~Type_of_cost,data=AustraliaHigh,FUN=sum)$cost_bil)
aggregate(cost_bil~Type_of_cost,data=AustraliaHigh,FUN=length)

# sector
aggregate(cost_bil~Impacted_sector,data=Australia,FUN="sum")
aggregate(cost_bil~Impacted_sector,data=Australia,FUN="length")

aggregate(cost_bil~Impacted_sector,data=AustraliaHigh,FUN="sum")
aggregate(cost_bil~Impacted_sector,data=AustraliaHigh,FUN="sum")$cost_bil / (sum(aggregate(cost_bil~Impacted_sector,data=AustraliaHigh,FUN="sum")$cost_bil))
aggregate(cost_bil~Impacted_sector,data=AustraliaHigh,FUN="length")


# by state (all)
# AUS = MULTIPLE STATES IN AUSTRALIA
# AUSTERR = AUSTRALIAN TERRITORY (NOT IN A STATE PER SE)
# ACT = AUSTRALIAN CAPITAL TERRITORY
# NSW = NEW SOUTH WALES
# QLD = QUEENSLAND
# SA = SOUTH AUSTRALIA
# TAS = TASMANIA
# VIC = VICTORIA
# WA = WESTERN AUSTRALIA

Australia$State <- revalue(Australia$State, c("Australia" = "AUS",
                                              "Australian Capital Territory" = "ACT",
                                              "Barrington Tops National Park"="NSW",
                                              "Barrow Island"="WA",
                                              "Berringan" = "NSW",
                                              "Booderee National Park"="NSW",
                                              "Brisbane"="QLD",
                                              "Brisbane/Ipswich "="QLD",
                                              "Bulman community"="NT",
                                              "Bush Heritage rangelands"="WA",
                                              "Cairns area"="QLD",
                                              "Calga Springs Sanctuary/Tidbinbilla Nature Reserve/Woodlands Historic Park and Hamilton Community Parklands/Cleland Wildlife Park"="AUS",
                                              "Canning River"="WA",
                                              "Canning River, Western Australia"="WA",
                                              "Christmas Island"="AUSTERR",
                                              "Christmas Island National Park"="AUSTERR",
                                              "Cobourg Peninsula (Garig Gunak Barlu National Park)"="NT",
                                              "Cocos (Keeling) Islands" = "AUSTERR",
                                              "Collinsvale, Qld" = "QLD",
                                              "Crown Land"="AUS",
                                              "Darling Downs/Queensland"="QLD",
                                              "Darwin area" = "NT",
                                              "Dirk Hartog Island\nArea (62 790 Ha)"="WA",
                                              "Diverse / Unspecified"="AUS",
                                              "El Arish and Japoon National Park, North Queensland"="QLD",
                                              "Faure Island"="WA",
                                              "Fitz-Sterling region (Western Australia)"="WA",
                                              "Fitz-Sterling State (Western Australia)"="WA",
                                              "Hermite"="WA",
                                              "Hunter Valley"="NSW",
                                              "Hypothetical Commercial Cattle Farm In Central Queensland Brigalow Belt"="QLD",
                                              "Illawarra"="NSW",
                                              "Kakadu National Park"="NT",
                                              "Kakadu National Park (Jim Jim district)"="NT",
                                              "Kakadu National Park (Mary River district)"="NT",
                                              "Kakadu National Park (South Alligator district)"="NT",
                                              "Kalka"="SA",
                                              "Kimberly"="WA",
                                              "Kosciusko National Park "="NSW",
                                              "Lake Eyre Basin"="WA",
                                              "Lake Eyre Bassin"="WA",
                                              "Little Desert Lodge and Malleefowl Sanctuary/Scotia Sanctuary/Venus Bay Conservation Park/Peron Peninsula/Paruna Sanctuary/Ellenbrook Nature Reserve/Twin Swamps Nature Reserve"="AUS",
                                              "Lord Howe Island"="AUSTERR",
                                              "Lorna Glen (Matuwa Conservation Park, Western Australia)"="WA",
                                              "Mackay region" = "QLD",
                                              "Macquarie Island"="AUSTERR",
                                              "Maroochy shire council" = "QLD",
                                              "Monto and Wide Bay burnett region" = "QLD",
                                              "Mulligan Flat"="ACT",
                                              "Muloorina Station"="SA",
                                              "Murray Valley" = "AUS",
                                              "New South Wales "="NSW",
                                              "New South Wales"="NSW",
                                              "New South Wales/Victoria" = "AUS",
                                              "North Queensland (11 banana farms/19 sugar cane farms)"="QLD",
                                              "Northern Australia"="AUS",
                                              "Northern Territory"="NT",
                                              "Northwern Queensland"="QLD",
                                              "Northwestern Australia"="WA",
                                              "NSW"="NSW",
                                              "NSW/QLD"="AUS",
                                              "NSW/Victoria/Tas" = "AUS",
                                              "NSW national parks"="NSW",
                                              "NT"="NT",
                                              "Paroo Shire"="QLD",
                                              "Phillip Island"="VIC", #yes mine... :-)
                                              "Pilbara Islands"="WA",
                                              "Private forests"="AUS",
                                              "Qld"="QLD",
                                              "Qld/NSW" = "AUS",
                                              "Queensland (Bulloo Downs Station)"="QLD",
                                              "Queensland (Northern Hairy-nosed wombat population range)"="QLD",
                                              "Queensland And New South Wales"="AUS",
                                              "Queensland, Peter Faust Dam"="QLD",
                                              "Queensland, Victoria And New South Wales"="AUS",
                                              "Queensland"="QLD",
                                              "Queensland'"="QLD",
                                              "Queensland/New South Wales"="AUS",
                                              "Queensland/NSW"="AUS",
                                              "Stateal Nt"="NT",
                                              "Regional, Nt"="NT",
                                              "Roadside, Midlands Highway"="TAS",
                                              "SA"="SA",
                                              "SA/NSW/WA/QLD" = "AUS",
                                              "Sandy Lacepede"="WA",
                                              "Scotia Sanctuary "="NSW",
                                              "South Australia"="SA",
                                              "South Australia (Pastoral State)"="SA",
                                              "Southern Australia"="AUS",
                                              "South Australia (Pastoral region)"="SA",
                                              "Southern half of Western Australia"="WA",
                                              "Sydney - Liverpool city council" = "NSW",
                                              "TA"="TAS",
                                              "Tas"="TAS",
                                              "Victoria/South Australia" = "AUS",
                                              "Tasmania"="TAS",
                                              "Tiwi Islands"="NT",
                                              "Tjukurla"="WA",
                                              "Vic"="VIC",
                                              "Vic/NSW" = "AUS",
                                              "Victoria"="VIC",
                                              "WA"="WA",
                                              "WA, SA, NT"="AUS",
                                              "Wastern Australia" = "WA",
                                              "Warakurna"="WA",
                                              "Wardang Island/Currawinya National Park/Royal Botanic Gardens Cranbourne/Heirisson Prong/Watarrka National Park/Yaraandoo Environmental Centre/N.S.W Roads and Traffic Authority"="AUS",
                                              "Western Australia "="WA",
                                              "Western Australia"="WA",
                                              "Western Division of New South Wales (about 877 000ha of national parks and reserves)"="NSW",
                                              "Wyperfeld National Park"="VIC",
                                              "Yookamurra reserve"="SA",
                                              "Yookamurra Sanctuary/Living Desert Wildlife Park/Karakamia"="AUS"))

# check
table(Australia$State)

State <- aggregate(cost_bil~State,data=Australia,FUN="sum")
State[order(State$cost_bil),]
State[order(State$cost_bil),]$cost_bil / (sum(State$cost_bil) - max(State$cost_bil))
(sum(State$cost_bil) - max(State$cost_bil))
aggregate(cost_bil~State,data=Australia,FUN="length")

lState <- aggregate(cost_bil~State,data=Australia,FUN="length")
colnames(lState) <- c("State","len")
State.mrg <- merge(State, lState, by='State') 

## relationships with # entries
all.states <- log10(State.mrg[,c(2:3)])
plot(all.states[,2],all.states[,1],pch=19)
lfit.all <- lm(all.states[,1] ~ all.states[,2])
summary(lfit.all)
abline(lfit.all,col="red",lty=2)
linreg.ER(all.states[,2],all.states[,1])
coef(lfit.all)
entries.vec <- seq(min(all.states$len),  max(all.states$len), by=0.01)
costs.pred <- 10^(coef(lfit.all)[1]+coef(lfit.all)[2]*entries.vec)
all.out <- data.frame(10^entries.vec,costs.pred)
colnames(all.out) <- c("entries","predcost")

## by area
StateA <- c(2358,7688287,801150,2527013,1729742,1347791,984321,227444,68401)
StateLab <- c("ACT","AUS","NSW","WA","QLD","NT","SA","VIC","TAS")
StateAdat <- data.frame(StateLab,StateA)
colnames(StateAdat) <- c("State","area")
StateAdatAllmrg <- merge(State.mrg,StateAdat, by="State")
StateAdatAllmrg$entriesA <- StateAdatAllmrg$len/StateAdatAllmrg$area
StateAdatAllmrg

# relationships
plot(log10(StateAdatAllmrg$area),log10(StateAdatAllmrg$cost_bil), pch=19)
lAfit <- lm(log10(cost_bil) ~ log10(area), data=StateAdatAllmrg)
linreg.ER(log10(StateAdatAllmrg$area),log10(StateAdatAllmrg$cost_bil))
abline(lAfit,lty=2,col="grey")
area.vec <- seq(min(StateAdatAllmrg$area),  max(StateAdatAllmrg$area), by=100)
costApred <- 10^(coef(lAfit)[1]+coef(lAfit)[2]*log10(area.vec))
Acost.dat <- data.frame(area.vec,costApred)
colnames(Acost.dat) <- c("area","predcost")

plot(log10(StateAdatAllmrg$area),log10(StateAdatAllmrg$len), pch=19)
llfit <- lm(log10(len) ~ log10(area), data=StateAdatAllmrg)
linreg.ER(log10(StateAdatAllmrg$area),log10(StateAdatAllmrg$len))
abline(llfit,lty=2,col="grey")
lenApred <- 10^(coef(llfit)[1]+coef(llfit)[2]*log10(area.vec))
Al.dat <- data.frame(area.vec,lenApred)
colnames(Al.dat) <- c("area","predlen")

# by state (reliable, observed)
StateHigh <- aggregate(cost_bil~State,data=AustraliaHigh,FUN="sum")
StateHigh[order(StateHigh$cost_bil),]
sum(StateHigh[order(StateHigh$cost_bil),]$cost_bil) - max(StateHigh[order(StateHigh$cost_bil),]$cost_bil)
aggregate(cost_bil~State,data=AustraliaHigh,FUN="length")

lStateHigh <- aggregate(cost_bil~State,data=AustraliaHigh,FUN="length")
colnames(lStateHigh) <- c("State","len")
StateHigh.mrg <- merge(StateHigh, lStateHigh, by='State') 

StateAdatHighmrg <- merge(StateHigh.mrg,StateAdat, by="State")
StateAdatHighmrg$entriesA <- StateAdatHighmrg$len/StateAdatHighmrg$area

# relationships
plot(log10(StateAdatHighmrg$area),log10(StateAdatHighmrg$cost_bil), pch=19)
lAfitHigh <- lm(log10(cost_bil) ~ log10(area), data=StateAdatHighmrg)
linreg.ER(log10(StateAdatHighmrg$area),log10(StateAdatHighmrg$cost_bil))
abline(lAfitHigh,lty=2,col="grey")
area.vec <- seq(min(StateAdatHighmrg$area),  max(StateAdatHighmrg$area), by=100)
costApredHigh <- 10^(coef(lAfitHigh)[1]+coef(lAfitHigh)[2]*log10(area.vec))
AcostHigh.dat <- data.frame(area.vec,costApredHigh)
colnames(AcostHigh.dat) <- c("area","predcost")

plot(log10(StateAdatHighmrg$area),log10(StateAdatHighmrg$len), pch=19)
llfitHigh <- lm(log10(len) ~ log10(area), data=StateAdatHighmrg)
linreg.ER(log10(StateAdatHighmrg$area),log10(StateAdatHighmrg$len))
abline(llfitHigh,lty=2,col="grey")
lenApredHigh <- 10^(coef(llfitHigh)[1]+coef(llfitHigh)[2]*log10(area.vec))
AlHigh.dat <- data.frame(area.vec,lenApredHigh)
colnames(AlHigh.dat) <- c("area","predlen")

reliable.states <- log10(StateHigh.mrg[,c(2:3)])
plot(reliable.states[,2],reliable.states[,1],pch=19)
lfit.reliable <- lm(reliable.states[,1] ~ reliable.states[,2])
summary(lfit.reliable)
abline(lfit.reliable,col="red",lty=2)
linreg.ER(reliable.states[,2],reliable.states[,1])
coef(lfit.reliable)
entries.vec <- seq(min(reliable.states$len),  max(reliable.states$len),  by=0.01)
costs.pred <- 10^(coef(lfit.reliable)[1]+coef(lfit.reliable)[2]*entries.vec)
reliable.out <- data.frame(10^entries.vec,costs.pred)
colnames(reliable.out) <- c("entries","predcost")

# most important species by state
# ACT
ACTHigh <- AustraliaHigh[AustraliaHigh$State %in% c("ACT"),]
ACT.spp <- aggregate(cost_bil~Species,data=ACTHigh,FUN="sum")
ACT.spp[order(ACT.spp$cost_bil),]
100 * ACT.spp[order(ACT.spp$cost_bil),]$cost_bil / sum(ACT.spp[order(ACT.spp$cost_bil),]$cost_bil)

# NSW
NSWHigh <- AustraliaHigh[AustraliaHigh$State %in% c("NSW"),]
NSW.spp <- aggregate(cost_bil~Species,data=NSWHigh,FUN="sum")
NSW.spp[order(NSW.spp$cost_bil),]
sum((NSW.spp[order(NSW.spp$cost_bil),][1:44,])$cost_bil)
100 * NSW.spp[order(NSW.spp$cost_bil),]$cost_bil / sum(NSW.spp[order(NSW.spp$cost_bil),]$cost_bil)

# NT
NTHigh <- AustraliaHigh[AustraliaHigh$State %in% c("NT"),]
NT.spp <- aggregate(cost_bil~Species,data=NTHigh,FUN="sum")
NT.spp[order(NT.spp$cost_bil),]
100 * NT.spp[order(NT.spp$cost_bil),]$cost_bil / sum(NT.spp[order(NT.spp$cost_bil),]$cost_bil)

# QLD
QLDHigh <- AustraliaHigh[AustraliaHigh$State %in% c("QLD"),]
QLD.spp <- aggregate(cost_bil~Species,data=QLDHigh,FUN="sum")
QLD.spp[order(QLD.spp$cost_bil),]
QLD.spp[order(QLD.spp$cost_bil),][1:33,]
sum(QLD.spp[order(QLD.spp$cost_bil),][1:33,]$cost_bil)

100 * QLD.spp[order(QLD.spp$cost_bil),]$cost_bil / sum(QLD.spp[order(QLD.spp$cost_bil),]$cost_bil)

# SA
SAHigh <- AustraliaHigh[AustraliaHigh$State %in% c("SA"),]
SA.spp <- aggregate(cost_bil~Species,data=SAHigh,FUN="sum")
SA.spp[order(SA.spp$cost_bil),]
sum(SA.spp[order(SA.spp$cost_bil),]$cost_bil[1:38])
100 * SA.spp[order(SA.spp$cost_bil),]$cost_bil / sum(SA.spp[order(SA.spp$cost_bil),]$cost_bil)

# TAS
TASHigh <- AustraliaHigh[AustraliaHigh$State %in% c("TAS"),]
TAS.spp <- aggregate(cost_bil~Species,data=TASHigh,FUN="sum")
TAS.spp[order(TAS.spp$cost_bil),]
sum(TAS.spp[order(TAS.spp$cost_bil),][1:19,]$cost_bil)
100 * TAS.spp[order(TAS.spp$cost_bil),]$cost_bil / sum(TAS.spp[order(TAS.spp$cost_bil),]$cost_bil)

# VIC
VICHigh <- AustraliaHigh[AustraliaHigh$State %in% c("VIC"),]
VIC.spp <- aggregate(cost_bil~Species,data=VICHigh,FUN="sum")
VIC.spp[order(VIC.spp$cost_bil),]
sum(VIC.spp[order(VIC.spp$cost_bil),][1:29,]$cost_bil)
100 * VIC.spp[order(VIC.spp$cost_bil),]$cost_bil / sum(VIC.spp[order(VIC.spp$cost_bil),]$cost_bil)

# WA
WAHigh <- AustraliaHigh[AustraliaHigh$State %in% c("WA"),]
WA.spp <- aggregate(cost_bil~Species,data=WAHigh,FUN="sum")
WA.spp[order(WA.spp$cost_bil),]
sum(WA.spp[order(WA.spp$cost_bil),][1:42,]$cost_bil)
100 * WA.spp[order(WA.spp$cost_bil),]$cost_bil / sum(WA.spp[order(WA.spp$cost_bil),]$cost_bil)


## Kingdom
Australia$Kingdom <- revalue(Australia$Kingdom, c("Animalia/Plantae" = "Diverse/Unspecified"))
Australia$Kingdom <- revalue(Australia$Kingdom, c("Orthornavirae" = "Viruses"))

# check
table(Australia$Kingdom)

Kingdom <- aggregate(cost_bil~Kingdom,data=Australia,FUN="sum")
Kingdom[order(Kingdom$cost_bil),]
100 * Kingdom[order(Kingdom$cost_bil),]$cost_bil / sum(Kingdom[order(Kingdom$cost_bil),]$cost_bil)
aggregate(cost_bil~Kingdom,data=Australia,FUN="length")

KingdomHigh <- aggregate(cost_bil~Kingdom,data=AustraliaHigh,FUN="sum")
KingdomHigh[order(KingdomHigh$cost_bil),]
100 * KingdomHigh[order(KingdomHigh$cost_bil),]$cost_bil / sum(KingdomHigh[order(KingdomHigh$cost_bil),]$cost_bil)
aggregate(cost_bil~Kingdom,data=AustraliaHigh,FUN="length")

KingdomLow <- aggregate(cost_bil~Kingdom,data=AustraliaLow,FUN="sum")
KingdomLow[order(KingdomLow$cost_bil),]
aggregate(cost_bil~Kingdom,data=AustraliaLow,FUN="length")


# Class
Australia$Class <- revalue(Australia$Class, c("Amphibia/Arthropoda/Mammalia/Diverse" = "Diverse/Unspecified",
                                              "Dicotyledonae" = "Eudicots",
                                              "Magnoliopsida" = "Eudicots",
                                              "Monocotyledonae" = "Monocots",
                                              "Liliopsida" = "Monocots"))

Australia$Class <- factor(Australia$Class)

# check
table(Australia$Class)

ClassAll <- aggregate(cost_bil~Class,data=Australia,FUN="sum")

ClassAll[order(ClassAll$cost_bil),]
aggregate(cost_bil~Class,data=Australia,FUN="length")
table(AustraliaHigh$Class)

AustraliaHigh$Class <- revalue(AustraliaHigh$Class, c("Liliopsida" = "Monocots",
                                              "Monocotyledonae"="Monocots",
                                              "Magnoliopsida"="Eudicots",
                                              "Dicotyledonae"="Eudicots"))

ClassHigh <- aggregate(cost_bil~Class,data=AustraliaHigh,FUN="sum")
ClassHigh[order(ClassHigh$cost_bil),]
aggregate(cost_bil~Class,data=AustraliaHigh,FUN="length")


## by specific taxonomic group
# plants
PlantsAll <- Australia[Australia$Kingdom %in% c("Plantae"),]
PlantsHigh <- AustraliaHigh[AustraliaHigh$Kingdom %in% c("Plantae"),]

PlantsAll$Species <- revalue(PlantsAll$Species, c("Kalanchoe delagoensis/Bryophyllum houghtonii/Bryophyllum pinnatum" = "diverse",
                                                      "Parkinsonia aculeata//Vachellia farnesiana" = "diverse",
                                                      "Carduus pycnocephalus and Carduus tenuiflorus" = "Carduus spp.",
                                                      "Andropogon gayanus/Cenchrus polystachios/Hymenachne amplexicaulis"="diverse",
                                                      "Acaciella angustissima/Aeschynomene paniculata/Aeschynomene brasiliana/Indigofera schimperi"="diverse",
                                                      "Acacia nilotica/Parkinsonia aculeata/Jatropha gossypiifolia/Ziziphus mauritiana"="diverse",
                                                      "Sporobolus pyramidalis/Sporobolus fertilis"="Sporobolus spp.",
                                                      "Proposis pallida/Proposis gladulosa/ Proposis juliflora/ Proposis velutina"="Prosopis spp.",
                                                      "Echium plantagineum"="Echium spp.",
                                                      "Parkinsonia aculeata/Cryptostegia grandiflora/Calotropis procera/Ziziphus mauritiana"="diverse",
                                                      "Echium plantagineum/Echium vulgare/Echium/italicum/Echium simplex"="Echium spp.",
                                                      "Sporobolus pyramidalis"="Sporobolus spp.",
                                                      "Emex australis and Emex spinosa" = "Emex spp.",
                                                      "Emex australis" = "Emex spp.",
                                                       "Pelargonium\xcaalchemilloides" = "Pelargonium spp.",
                                                      "Sida acuta and Sida rhombifolia" = "Sida spp.",
                                                      "Amaranthus\xcamitchellii" = "Amaranthus spp.",
                                                      "Sporobolus anglicus"="Sporobolus spp.",
                                                      "Salvinia molesta"="Salvinia spp.",
                                                      "Diverse/Unspecified"="diverse"))

# check
table(PlantsAll$Species)

PlantsHigh$Species <- revalue(PlantsHigh$Species, c("Kalanchoe delagoensis/Bryophyllum houghtonii/Bryophyllum pinnatum" = "diverse",
                                                  "Parkinsonia aculeata//Vachellia farnesiana" = "diverse",
                                                  "Amaranthus\xcamitchellii" = "Amaranthus spp.",
                                                  "Andropogon gayanus/Cenchrus polystachios/Hymenachne amplexicaulis"="diverse",
                                                  "Acaciella angustissima/Aeschynomene paniculata/Aeschynomene brasiliana/Indigofera schimperi"="diverse",
                                                  "Acacia nilotica/Parkinsonia aculeata/Jatropha gossypiifolia/Ziziphus mauritiana"="diverse",
                                                  "Carduus pycnocephalus and Carduus tenuiflorus" = "Carduus spp.",
                                                  "Emex australis and Emex spinosa" = "Emex spp.",
                                                  "Emex australis" = "Emex spp.",
                                                  "Sporobolus pyramidalis/Sporobolus fertilis"="Sporobolus spp.",
                                                  "Proposis pallida/Proposis gladulosa/ Proposis juliflora/ Proposis velutina"="Prosopis spp.",
                                                  "Echium plantagineum"="Echium spp.",
                                                  "Parkinsonia aculeata/Cryptostegia grandiflora/Calotropis procera/Ziziphus mauritiana"="diverse",
                                                  "Echium plantagineum/Echium vulgare/Echium/italicum/Echium simplex"="Echium spp.",
                                                  "Sporobolus pyramidalis"="Sporobolus spp.",
                                                  "Pelargonium\xcaalchemilloides" = "Pelargonium spp.",
                                                  "Sporobolus anglicus"="Sporobolus spp.",
                                                  "Sida acuta and Sida rhombifolia" = "Sida spp.",
                                                  "Salvinia molesta"="Salvinia spp.",
                                                  "Diverse/Unspecified"="diverse"))
# check
table(PlantsHigh$Species)

PlantsAll.spp <- aggregate(cost_bil~Species, data=PlantsAll, FUN="sum")
PlantsAll.spp.sort <- PlantsAll.spp[order(PlantsAll.spp$cost_bil),]
PlantsAll.spp.sort
sum(PlantsAll.spp.sort$cost_bil)
PlantsAll.spp.sort[1:(length(PlantsAll.spp.sort$Species)-21),]
sum(PlantsAll.spp.sort[1:(length(PlantsAll.spp.sort$Species)-21),]$cost_bil)

sum(PlantsAll.spp[order(PlantsAll.spp$cost_bil),]$cost_bil)
length(PlantsAll.spp[order(PlantsAll.spp$cost_bil),]$cost_bil)


PlantsHigh.spp <- aggregate(cost_bil~Species, data=PlantsHigh, FUN="sum")
PlantsHigh.spp
PlantsHigh.spp.sort <- PlantsAll.spp[order(PlantsAll.spp$cost_bil),]
PlantsHigh.spp.sort
PlantsHigh.spp.sort[1:(length(PlantsHigh.spp.sort$Species)-21),]
sum(PlantsHigh.spp.sort[1:(length(PlantsHigh.spp.sort$Species)-21),]$cost_bil)
sum(PlantsHigh.spp[order(PlantsHigh.spp$cost_bil),]$cost_bil)
length(PlantsHigh.spp[order(PlantsHigh.spp$cost_bil),]$cost_bil)

sum(PlantsHigh.spp[order(PlantsHigh.spp$cost_bil),]$cost_bil)
data.frame(PlantsHigh.spp.sort$Species, 100*(PlantsHigh.spp.sort$cost_bil)/(sum(PlantsHigh.spp[order(PlantsHigh.spp$cost_bil),]$cost_bil)))

length(PlantsHigh.spp[order(PlantsHigh.spp$cost_bil),]$cost_bil)

PlantsHigh.type <- aggregate(cost_bil~Type_of_cost, data=PlantsHigh,FUN="sum")

# Magnoliopsida
Magnols <- Australia[Australia$Class %in% c("Magnoliopsida"),]
Magnols.spp <- aggregate(cost_bil~Species,data=Magnols,FUN="sum")
Magnols.spp[order(Magnols.spp$cost_bil),]
100 * Magnols.spp[order(Magnols.spp$cost_bil),]$cost_bil / sum(Magnols.spp[order(Magnols.spp$cost_bil),]$cost_bil)

# Liliopsida
Lilis <- Australia[Australia$Class %in% c("Liliopsida"),]
Lilis.spp <- aggregate(cost_bil~Species,data=Lilis,FUN="sum")
Lilis.spp[order(Lilis.spp$cost_bil),]
100 * Lilis.spp[order(Lilis.spp$cost_bil),]$cost_bil / sum(Lilis.spp[order(Lilis.spp$cost_bil),]$cost_bil)

# Eudicots
Eudicots <- Australia[Australia$Class %in% c("Eudicots"),]
Eudicots.spp <- aggregate(cost_bil~Species,data=Eudicots,FUN="sum")
Eudicots.spp[order(Eudicots.spp$cost_bil),]
100 * Eudicots.spp[order(Eudicots.spp$cost_bil),]$cost_bil / sum(Eudicots.spp[order(Eudicots.spp$cost_bil),]$cost_bil)

# Monocots
Monocots <- Australia[Australia$Class %in% c("Monocots"),]
table(Monocots$Class)
Monocots.spp <- aggregate(cost_bil~Species,data=Monocots,FUN="sum")
Monocots.spp[order(Monocots.spp$cost_bil),]
100 * Monocots.spp[order(Monocots.spp$cost_bil),]$cost_bil / sum(Monocots.spp[order(Monocots.spp$cost_bil),]$cost_bil)

# Ulvophyceae
Ulvophytes <- Australia[Australia$Class %in% c("Ulvophyceae"),]
Ulvophytes.spp <- aggregate(cost_bil~Species,data=Ulvophytes,FUN="sum")
Ulvophytes.spp[order(Ulvophytes.spp$cost_bil),]

# brown algae
BrownAlgae <- Australia[Australia$Class %in% c("Phaeophyceae"),]
BrownAlgae.spp <- aggregate(cost_bil~Species,data=BrownAlgae,FUN="sum")
BrownAlgae.spp[order(BrownAlgae.spp$cost_bil),]

# ferns
Ferns <- Australia[Australia$Class %in% c("Polypodiopsida"),]
Ferns.spp <- aggregate(cost_bil~Species,data=Ferns,FUN="sum")
Ferns.spp[order(Ferns.spp$cost_bil),]

# fungi
# Dothideomycetes
Dfungi <- Australia[Australia$Class %in% c("Dothideomycetes"),]
Dfungi.spp <- aggregate(cost_bil~Species,data=Dfungi,FUN="sum")
Dfungi.spp[order(Dfungi.spp$cost_bil),]

# Sordariomycetes
Sorda <- Australia[Australia$Class %in% c("Sordariomycetes"),]
Sorda.spp <- aggregate(cost_bil~Species,data=Sorda,FUN="sum")
Sorda.spp[order(Sorda.spp$cost_bil),]

# Agaricomycetes
Agarics <- Australia[Australia$Class %in% c("Agaricomycetes"),]
Agarics.spp <- aggregate(cost_bil~Species,data=Agarics,FUN="sum")
Agarics.spp[order(Agarics.spp$cost_bil),]

# Pucciniomycetes
Puccin <- Australia[Australia$Class %in% c("Pucciniomycetes"),]
Puccin.spp <- aggregate(cost_bil~Species,data=Puccin,FUN="sum")
Puccin.spp[order(Puccin.spp$cost_bil),]

# Leotiomycetes
Leotio <- Australia[Australia$Class %in% c("Leotiomycetes"),]
Leotio.spp <- aggregate(cost_bil~Species,data=Leotio,FUN="sum")
Leotio.spp[order(Leotio.spp$cost_bil),]

# Exobasidiomycetes
Exobasids <- Australia[Australia$Class %in% c("Exobasidiomycetes"),]
Exobasids.spp <- aggregate(cost_bil~Species,data=Exobasids,FUN="sum")
Exobasids.spp[order(Exobasids.spp$cost_bil),]

# bacteria
GPB <- Australia[Australia$Class %in% c("Gammaproteobacteria"),]
GPB.spp <- aggregate(cost_bil~Species,data=GPB,FUN="sum")
GPB.spp[order(GPB.spp$cost_bil),]

# nematodes
Nematodes <- Australia[Australia$Class %in% c("Secernentea"),]
Nematodes.spp <- aggregate(cost_bil~Species,data=Nematodes,FUN="sum")
Nematodes.spp[order(Nematodes.spp$cost_bil),]

# gastropods
Gastropods <- Australia[Australia$Class %in% c("Gastropoda"),]
Gastropods.spp <- aggregate(cost_bil~Species,data=Gastropods,FUN="sum")
Gastropods.spp[order(Gastropods.spp$cost_bil),]

# insects
Insects <- Australia[Australia$Class %in% c("Insecta"),]
Insects$Species <- revalue(Insects$Species, c("Anoplolepis gracilipes"="Anoplolepis gracilipes",
                                              "Anoploepsis gracilipes"="Anoplolepis gracilipes"))

Insects.spp <- aggregate(cost_bil~Species,data=Insects,FUN="sum")
Insects.spp[order(Insects.spp$cost_bil),]

InsectsObs <- Insects[Insects$Implementation %in% c("Observed"),]
InsectsHigh <- InsectsObs[InsectsObs$Method_reliability %in% c("High"),]

InsectsHigh.spp <- aggregate(cost_bil~Species,data=InsectsHigh,FUN="sum")
InsectsHigh.spp[order(InsectsHigh.spp$cost_bil),]
100 * InsectsHigh.spp[order(InsectsHigh.spp$cost_bil),]$cost_bil / sum(InsectsHigh.spp[order(InsectsHigh.spp$cost_bil),]$cost_bil)
100 * InsectsHigh.spp[order(InsectsHigh.spp$cost_bil),]$cost_bil[-length(InsectsHigh.spp$Species)] / sum(InsectsHigh.spp[order(InsectsHigh.spp$cost_bil),]$cost_bil[-length(InsectsHigh.spp$Species)])
sum(100 * InsectsHigh.spp[order(InsectsHigh.spp$cost_bil),]$cost_bil[-length(InsectsHigh.spp$Species)] / sum(InsectsHigh.spp[order(InsectsHigh.spp$cost_bil),]$cost_bil[-length(InsectsHigh.spp$Species)]))

# arachnids
Arachnids <- Australia[Australia$Class %in% c("Arachnida"),]
Arachnids.spp <- aggregate(cost_bil~Species,data=Arachnids,FUN="sum")
Arachnids.spp[order(Arachnids.spp$cost_bil),]

# amphibians
Amphibians <- Australia[Australia$Class %in% c("Amphibia"),]
Amphibians.spp <- aggregate(cost_bil~Species,data=Amphibians,FUN="sum")
Amphibians.spp[order(Amphibians.spp$cost_bil),]
AmphibiansHigh <- AustraliaHigh[AustraliaHigh$Class %in% c("Amphibia"),]

# birds
Birds <- Australia[Australia$Class %in% c("Aves"),]
aggregate(cost_bil~Species,data=Birds,FUN="sum")

BirdsHigh <- AustraliaHigh[AustraliaHigh$Class %in% c("Aves"),]
aggregate(cost_bil~Species,data=BirdsHigh,FUN="sum")

# mammals
# mammals
Mammals <- Australia[Australia$Class %in% c("Mammalia"),]
Mammals$Common_name <- revalue(Mammals$Common_name, c("Black rat" = "rodent",
                                                      "Buffalo" = "buffalo",
                                                      "Camel"="camel",
                                                      "Cat/Red fox"="cat/fox",
                                                      "Cat"="cat",
                                                      "Cat/Red fox/Rabbit"="cat/fox/rabbit",
                                                      "Deer"="deer",
                                                      "Deers"="deer",
                                                      "Dingo"="dingo",
                                                      "Dingo/Dog"="dingo",
                                                      "Dingo/Dog/European rabbits/Feral pig"="dingo/rabbit/pig",
                                                      "Dingo/Dog/Red fox"="dingo/fox",
                                                      "Dingo/feral dog"="dingo",
                                                      "Dingo/Wild dogs"="dingo",
                                                      "Diverse (may include non invasive pests)"="diverse",
                                                      "Dog"="dingo",
                                                      "Donkey"="equid",
                                                      "Donkeys"="equid",
                                                      "Dromedary camel"="camel",
                                                      "European (common) rabbit"="rabbit",
                                                      "European rabbit/Rodents"="rabbit/rodent",
                                                      "Feral Cat"="cat",
                                                      "Feral cats"="cat",
                                                      "Feral cats/Dingoes"="cat/dingo",
                                                      "Feral goat"="goat",
                                                      "Feral goats"="goat",
                                                      "Feral pig"="pig",
                                                      "Feral pigs"="pig",
                                                      "Foxes"="fox",
                                                      "Goat"="goat",
                                                      "Horse"="horse",
                                                      "Horses "="horse",
                                                      "horse/donkey"="equid",
                                                      "Horses/donkeys"="equid",
                                                      "House mouse"="rodent",
                                                      "House mouse/European rabbit/Black rat"="rabbit/rodent",
                                                      "Koala"="koala",
                                                      "Macropod"="kangaroo",
                                                      "Mammalian pests"="diverse",
                                                      "Mammals"="diverse",
                                                      "Mice"="rodent",
                                                      "Mice and rats"="rodent",
                                                      "Pig"="pig",
                                                      "Predators"="cat/fox",
                                                      "Predators (cats and foxes)"="cat/fox",
                                                      "Rabbits"="rabbit",
                                                      "Rabbit"="rabbit",
                                                      "Rat"="rodent",
                                                      "Rat/Mouse"="rodent",
                                                      "Rats"="rodent",
                                                      "Red fox"="fox",
                                                      "Rodent"="rodent",
                                                      "Sheep"="sheep",
                                                      "Swamp buffalo"="buffalo",
                                                      "Water buffalo"="buffalo",
                                                      "Wild deers"="deer",
                                                      "Wild dogs"="dingo",
                                                      "Wild Horse"="horse",
                                                      "Wombat"="wombat"))

Mammals.spp <- aggregate(cost_bil~Common_name,data=Mammals,FUN="sum")
Mammals.spp[order(Mammals.spp$cost_bil),]
100 * Mammals.spp[order(Mammals.spp$cost_bil),]$cost_bil / sum(Mammals.spp[order(Mammals.spp$cost_bil),]$cost_bil)
sum(Mammals.spp[order(Mammals.spp$cost_bil),]$cost_bil)

Mammals.sppL <- aggregate(cost_bil~Species,data=Mammals,FUN="sum")
Mammals.sppL[order(Mammals.sppL$cost_bil),]


MammalsObs <- Mammals[Mammals$Implementation %in% c("Observed"),]
MammalsHigh <- MammalsObs[MammalsObs$Method_reliability %in% c("High"),]
MammalsHigh.spp <- aggregate(cost_bil~Common_name,data=MammalsHigh,FUN="sum")
MammalsHigh.spp[order(MammalsHigh.spp$cost_bil),]
100 * MammalsHigh.spp[order(MammalsHigh.spp$cost_bil),]$cost_bil / sum(MammalsHigh.spp[order(MammalsHigh.spp$cost_bil),]$cost_bil)
sum(MammalsHigh.spp[order(MammalsHigh.spp$cost_bil),]$cost_bil)

MammalsHigh.type <- aggregate(cost_bil~Type_of_cost, data=MammalsHigh,FUN="sum")



##################
## costs per year
##################

## all costs (NOTE: replace function 'calculateRawAvgCosts' with 'summarizeCosts' in latest version of R)
CostsPerYearRaw <- calculateRawAvgCosts(Australia, cost.column = "cost_bil",  in.millions = FALSE,  minimum.year = 1970, maximum.year = 2020)

CostsPerYearRaw$cost.per.year
CostsPerYearRaw$average.cost.per.period
CostsPerYearRaw$average.total.cost

## cumulative costs
cost.cum <- cumsum(aggregate(CostsPerYearRaw$cost.per.year[,2],by=list(CostsPerYearRaw$cost.per.year[,1]),sum, na.rm=T)[,2])
est.cum <- cumsum(aggregate(CostsPerYearRaw$cost.per.year[,3],by=list(CostsPerYearRaw$cost.per.year[,1]),sum, na.rm=T)[,2])
Raw.dat <- data.frame(est.cum[-c(1:2)],cost.cum[-c(1:2)])
colnames(Raw.dat) <- c("cumest","cumcost")

# back-transform
Raw.est.vec <- as.data.frame(seq(1,max(Raw.dat$cumest)))
colnames(Raw.est.vec) <- c("cumest")
pred.cost <- 10^(coef(lfit)[1] + coef(lfit)[2]*log10(Raw.est.vec$cumest))

all.cumcostall.out <- data.frame(Raw.dat$cumest, Raw.dat$cumcost)
colnames(all.cumcostall.out) <- c("cumest","cumcost")
all.cumcostpred.out <- data.frame(Raw.est.vec$cumest, pred.cost)
colnames(all.cumcostpred.out) <- c("cumest","predcumcost")


## reliable costs  (NOTE: replace function 'calculateRawAvgCosts' with 'summarizeCosts' in latest version of R)
CostsPerYearHigh <- calculateRawAvgCosts(AustraliaHigh, cost.column = "cost_bil",  in.millions = FALSE,  minimum.year = 1970, maximum.year = 2020)

CostsPerYearHigh$cost.per.year
CostsPerYearHigh$average.cost.per.period
CostsPerYearHigh$average.total.cost

## cumulative costs
cost.cum <- cumsum(aggregate(CostsPerYearHigh$cost.per.year[,2],by=list(CostsPerYearRaw$cost.per.year[,1]),sum, na.rm=T)[,2])
est.cum <- cumsum(aggregate(CostsPerYearHigh$cost.per.year[,3],by=list(CostsPerYearRaw$cost.per.year[,1]),sum, na.rm=T)[,2])
High.dat <- data.frame(est.cum[-c(1:2)],cost.cum[-c(1:2)])
colnames(High.dat) <- c("cumest","cumcost")

# back-transform
High.est.vec <- as.data.frame(seq(1,max(High.dat$cumest)))
colnames(High.est.vec) <- c("cumest")

pred.cost <- 10^(coef(lfit)[1] + coef(lfit)[2]*log10(High.est.vec$cumest))

high.cumcostraw.out <- data.frame(High.dat$cumest, High.dat$cumcost)
colnames(high.cumcostraw.out) <- c("cumest","cumcost")
high.cumcostpred.out <- data.frame(High.est.vec$cumest, pred.cost)
colnames(high.cumcostpred.out) <- c("cumest","predcumcost")


# by taxon (high only)
## Plants
CostsPerYearPlants <- calculateRawAvgCosts(PlantsAll, cost.column = "cost_bil",  in.millions = FALSE,  minimum.year = 1970, maximum.year = 2020)
CostsPerYearHighPlants <- calculateRawAvgCosts(PlantsHigh, cost.column = "cost_bil",  in.millions = FALSE,  minimum.year = 1970, maximum.year = 2020)
CostsPerYearHighPlants 
plot(CostsPerYearHighPlants)

# Mammals
CostsPerYearMammals <- calculateRawAvgCosts(Mammals, cost.column = "cost_bil",  in.millions = FALSE,  minimum.year = 1970, maximum.year = 2020)
CostsPerYearHighMammals <- calculateRawAvgCosts(MammalsHigh, cost.column = "cost_bil",  in.millions = FALSE,  minimum.year = 1970, maximum.year = 2020)
CostsPerYearHighMammals 
plot(CostsPerYearHighMammals)

# Insects
CostsPerYearInsects <- calculateRawAvgCosts(Insects, cost.column = "cost_bil",  in.millions = FALSE,  minimum.year = 1970, maximum.year = 2020)
CostsPerYearHighInsects <- calculateRawAvgCosts(InsectsHigh, cost.column = "cost_bil",  in.millions = FALSE,  minimum.year = 1970, maximum.year = 2020)
CostsPerYearHighInsects
plot(CostsPerYearHighInsects)



## time lag
# all costs
Australia$Publication_lag <- Australia$Publication_year - Australia$Impact_year

quantiles <- quantile(Australia$Publication_lag, probs = c(.25, .5, .75), na.rm=TRUE)
summary(quantiles)

# NOTE: replace function 'costTrendOverTime' with 'modelCosts' in latest version of R
cost.over.time <- costTrendOverTime(Australia,
                                    cost.column = "cost_bil",
                                    minimum.year = 1960,
                                    maximum.year = 2020,
                                    in.millions = FALSE,
                                    incomplete.year.threshold = 2020 - quantiles["75%"],
                                    gam.k = 8)

summarized.summary <- prettySummary(cost.over.time)
summarized.summary
cost.over.time$model.summary
cost.over.time$RMSE

RRl.slope <- cost.over.time$fitted.models$robust.linear$coefficients[2]
x.fold.incr <- 10^(RRl.slope * 10)

# OLS linear
AIC.OLSl <- AIC(cost.over.time$fitted.models$ols.linear)

# OLS quadratic
AIC.OLSq <- AIC(cost.over.time$fitted.models$ols.quadratic)

# GAM
AIC.GAM <- AIC(cost.over.time$fitted.models$gam)

AIC.vec <- c(AIC.OLSl, AIC.OLSq, AIC.GAM)
dAIC.vec <- delta.IC(AIC.vec)
AICwt.vec <- weight.IC(dAIC.vec)

fitted.costs <- cost.over.time$estimated.annual.costs
fitted.costs.OLSl <- subset(fitted.costs, model=="OLS regression" & Details == "Linear")[,c(2,4:6)]
fitted.costs.OLSq <- subset(fitted.costs, model=="OLS regression" & Details == "Quadratic")[,c(2,4:6)]
fitted.costs.GAM <- subset(fitted.costs, model=="GAM")[,c(2,4:6)]
fitted.costs.RRl <- subset(fitted.costs, model=="Robust regression" & Details == "Linear")[,c(2,4:6)]
fitted.costs.RRq <- subset(fitted.costs, model=="Robust regression" & Details == "Quadratic")[,c(2,4:6)]
fitted.costs.MARS <- subset(fitted.costs, model=="MARS")[,c(2,4:6)]

final.costs <- c(fitted.costs.OLSl$fit[length(fitted.costs.OLSl$fit)-3],
                 fitted.costs.OLSq$fit[length(fitted.costs.OLSq$fit)-3],
                 fitted.costs.RRl$fit[length(fitted.costs.RRl$fit)-3],
                 fitted.costs.RRq$fit[length(fitted.costs.RRq$fit)-3],
                 fitted.costs.GAM$fit[length(fitted.costs.GAM$fit)-3])

final.costs.lo <- c(fitted.costs.OLSl$lwr[length(fitted.costs.OLSl$lwr)-3],
                 fitted.costs.OLSq$lwr[length(fitted.costs.OLSq$lwr)-3],
                 fitted.costs.RRl$lwr[length(fitted.costs.RRl$lwr)-3],
                 fitted.costs.RRq$lwr[length(fitted.costs.RRq$lwr)-3],
                 fitted.costs.GAM$lwr[length(fitted.costs.GAM$lwr)-3])

final.costs.up <- c(fitted.costs.OLSl$upr[length(fitted.costs.OLSl$upr)-3],
                    fitted.costs.OLSq$upr[length(fitted.costs.OLSq$upr)-3],
                    fitted.costs.RRl$upr[length(fitted.costs.RRl$upr)-3],
                    fitted.costs.RRq$upr[length(fitted.costs.RRq$upr)-3],
                    fitted.costs.GAM$upr[length(fitted.costs.GAM$upr)-3])

wtAIC.final.costs <- AICwt.vec * final.costs[c(1,2,5)]
wtAIC.final.cost <- sum(wtAIC.final.costs)
wtAIC.final.cost

wtAIC.final.costs.lo <- AICwt.vec * final.costs.lo[c(1,2,5)]
wtAIC.final.cost.lo <- sum(wtAIC.final.costs.lo)
wtAIC.final.cost.lo

wtAIC.final.costs.up <- AICwt.vec * final.costs.up[c(1,2,5)]
wtAIC.final.cost.up <- sum(wtAIC.final.costs.up)
wtAIC.final.cost.up

dRMSE <- delta.IC(cost.over.time$RMSE[c(1:4,6),1])
RMSE.wt <- weight.IC(dRMSE)

wtRMSE.final.costs <- RMSE.wt * final.costs
wtRMSE.final.cost <- sum(wtRMSE.final.costs)
wtRMSE.final.cost

wtRMSE.final.costs.lo <- RMSE.wt * final.costs.lo
wtRMSE.final.cost.lo <- sum(wtRMSE.final.costs.lo)
wtRMSE.final.cost.lo

wtRMSE.final.costs.up <- RMSE.wt * final.costs.up
wtRMSE.final.cost.up <- sum(wtRMSE.final.costs.up)
wtRMSE.final.cost.up



## reliable costs
# lag
AustraliaHigh$Publication_lag <- AustraliaHigh$Publication_year - AustraliaHigh$Impact_year

quantiles <- quantile(AustraliaHigh$Publication_lag, probs = c(.25, .5, .75), na.rm=TRUE)
summary(quantiles)

cost.over.time <- costTrendOverTime(AustraliaHigh,
                                    cost.column = "cost_bil",
                                    minimum.year = 1960,
                                    maximum.year = 2020,
                                    in.millions = FALSE,
                                    incomplete.year.threshold = 2020 - quantiles["75%"],
                                    gam.k = 8)

cost.over.time
summarized.summary <- prettySummary(cost.over.time)
summarized.summary
cost.over.time$model.summary
cost.over.time$RMSE

RRl.slope <- cost.over.time$fitted.models$robust.linear$coefficients[2]
x.fold.incr <- 10^(RRl.slope * 10)

# OLS linear
AIC.OLSl <- AIC(cost.over.time$fitted.models$ols.linear)

# OLS quadratic
AIC.OLSq <- AIC(cost.over.time$fitted.models$ols.quadratic)

# GAM
AIC.GAM <- AIC(cost.over.time$fitted.models$gam)

AIC.vec <- c(AIC.OLSl, AIC.OLSq, AIC.GAM)
dAIC.vec <- delta.IC(AIC.vec)
AICwt.vec <- weight.IC(dAIC.vec)

fitted.costs <- cost.over.time$estimated.annual.costs
fitted.costs.OLSl <- subset(fitted.costs, model=="OLS regression" & Details == "Linear")[,c(2,4:6)]
fitted.costs.OLSq <- subset(fitted.costs, model=="OLS regression" & Details == "Quadratic")[,c(2,4:6)]
fitted.costs.GAM <- subset(fitted.costs, model=="GAM")[,c(2,4:6)]
fitted.costs.RRl <- subset(fitted.costs, model=="Robust regression" & Details == "Linear")[,c(2,4:6)]
fitted.costs.RRq <- subset(fitted.costs, model=="Robust regression" & Details == "Quadratic")[,c(2,4:6)]
fitted.costs.MARS <- subset(fitted.costs, model=="MARS")[,c(2,4:6)]

final.costs <- c(fitted.costs.OLSl$fit[length(fitted.costs.OLSl$fit)-3],
                 fitted.costs.OLSq$fit[length(fitted.costs.OLSq$fit)-3],
                 fitted.costs.RRl$fit[length(fitted.costs.RRl$fit)-3],
                 fitted.costs.RRq$fit[length(fitted.costs.RRq$fit)-3],
                 fitted.costs.GAM$fit[length(fitted.costs.GAM$fit)-3])

final.costs.lo <- c(fitted.costs.OLSl$lwr[length(fitted.costs.OLSl$lwr)-3],
                    fitted.costs.OLSq$lwr[length(fitted.costs.OLSq$lwr)-3],
                    fitted.costs.RRl$lwr[length(fitted.costs.RRl$lwr)-3],
                    fitted.costs.RRq$lwr[length(fitted.costs.RRq$lwr)-3],
                    fitted.costs.GAM$lwr[length(fitted.costs.GAM$lwr)-3])

final.costs.up <- c(fitted.costs.OLSl$upr[length(fitted.costs.OLSl$upr)-3],
                    fitted.costs.OLSq$upr[length(fitted.costs.OLSq$upr)-3],
                    fitted.costs.RRl$upr[length(fitted.costs.RRl$upr)-3],
                    fitted.costs.RRq$upr[length(fitted.costs.RRq$upr)-3],
                    fitted.costs.GAM$upr[length(fitted.costs.GAM$upr)-3])

wtAIC.final.costs <- AICwt.vec * final.costs[c(1,2,5)]
wtAIC.final.cost <- sum(wtAIC.final.costs)
wtAIC.final.cost

wtAIC.final.costs.lo <- AICwt.vec * final.costs.lo[c(1,2,5)]
wtAIC.final.cost.lo <- sum(wtAIC.final.costs.lo)
wtAIC.final.cost.lo

wtAIC.final.costs.up <- AICwt.vec * final.costs.up[c(1,2,5)]
wtAIC.final.cost.up <- sum(wtAIC.final.costs.up)
wtAIC.final.cost.up

dRMSE <- delta.IC(cost.over.time$RMSE[c(1:4,6),1])
RMSE.wt <- weight.IC(dRMSE)

wtRMSE.final.costs <- RMSE.wt * final.costs
wtRMSE.final.cost <- sum(wtRMSE.final.costs)
wtRMSE.final.cost

wtRMSE.final.costs.lo <- RMSE.wt * final.costs.lo
wtRMSE.final.cost.lo <- sum(wtRMSE.final.costs.lo)
wtRMSE.final.cost.lo

wtRMSE.final.costs.up <- RMSE.wt * final.costs.up
wtRMSE.final.cost.up <- sum(wtRMSE.final.costs.up)
wtRMSE.final.cost.up


