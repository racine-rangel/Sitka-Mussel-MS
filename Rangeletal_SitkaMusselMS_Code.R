################################################################
#Title: Rangel et al. Climate Manipulation Data Analysis
#Created by: R. Rangel
#Created: 05/2020
#Last edited: 04/2025
################################################################
##### Packages #####
install(plyr)
install(ggplot2)
install(RColorBrewer)
library(plyr)
library(ggplot2)
library(RColorBrewer)
library(lme4)
library(lmerTest)
library(nortest)
library(emmeans)
library(dplyr)
library(forcats)
library(reshape2)
library(pastecs) #Descriptive Stats
library(ggpubr) #Figures
library(forcats) #fct_relvel
library(pals) #for color in Community Comp figures
library(performance) #Check Model Residuals
library(seacarb) # carbonate chemistry
library(scales) # calculating % between strength values
#############################################
#Set working Directory
#############################################
setwd("/Users/racinerangel/Desktop/Climate Manipulations/Manuscript Data Files")
#######################################################

####################--------------------------------------------------------------
#Environmental Data
####################--------------------------------------------------------------


####################--------------------------------------------------------------
#Tide Pool Properties
####################--------------------------------------------------------------
Pools<- read.csv("TidePool_Propertities.csv")
str(Pools) #Check structure of data frame

Pools$Treatment<-factor(Pools$Treatment, levels = c("Unmanipulated", "CO2 Added", "Warmed", "CO2 Added + Warmed"))

TdHt_T <- Pools %>%
  group_by(Treatment) %>%
  summarize(mean = round(mean(Tide.Height..m.), 2),
            sd = round(sd(Tide.Height..m.), 2))

TdHt <- Pools %>%
  summarize(mean = round(mean(Tide.Height..m.), 2),
            sd = round(sd(Tide.Height..m.), 2))
# 2.56 m +/- 0.38

Vol <- Pools %>%
  summarize(mean = round(mean(Volume...pumped..L.), 2),
            sd = round(sd(Volume...pumped..L.), 2))
# 13.18 +/- 8.2 L

####################--------------------------------------------------------------
#pH Data All Pools - Hanna Data
####################--------------------------------------------------------------
pH_Meta<- read.csv("pH_Meta.csv")
str(pH_Meta) #Check structure of data frame

Day<-filter(pH_Meta, DayNight == "Day" )
Night<-filter(pH_Meta, DayNight == "Night")

Day$Month<-factor(Day$Month, levels = c("March", "July", "September"), labels = c("0", "3", "6"))
Day$Treatment<-factor(Day$Treatment, levels=c("Unmanipulated", "CO2 Added", "Warmed", "CO2 Added + Warmed"))

Night$Month<-factor(Night$Month, levels = c("March", "July", "September"), labels = c("0", "3", "6"))
Night$Treatment<-factor(Night$Treatment, levels=c("Unmanipulated", "CO2 Added", "Warmed", "CO2 Added + Warmed"))

####################--------------------------------------------------------------
#Descriptive stats by Month and Pool Treatment
####################--------------------------------------------------------------
Day%>% count(Treatment, Month)
####################--------------------------------------------------------------
#pH 
####################--------------------------------------------------------------
pH.36<-filter(pH_Meta, !Month =="March")

#CO2 pool difference between no CO2
pH.36 %>%
  group_by(CO2) %>%
  summarize(mean = mean(HannapH_Tot),
            sd = sd(HannapH_Tot),
            n = n(),
            se=sd/sqrt(n))
Day %>%
  group_by(Month, Treatment) %>%
  count(Treatment)

Night %>%
  group_by(Month, Treatment) %>%
  summarize(mean = mean(HannapH_Tot),
            sd = sd(HannapH_Tot))

####################--------------------------------------------------------------
#Temperature from Hanna
####################--------------------------------------------------------------
TempD <- Day %>%
  group_by(Month, Treatment) %>%
  summarize(mean = round(mean(Temp), 2),
            sd = round(sd(Temp), 2))
Temp_Day %>%
  group_by(Month, Treatment) %>%
  count(Treatment)

TempN <- Night %>%
  group_by(Month, Treatment) %>%
  summarize(mean = round(mean(Temp), 2),
            sd = round(sd(Temp), 2))

####################--------------------------------------------------------------
#Total Alkalinity
####################--------------------------------------------------------------
TADay<- Day %>%
  group_by(Month, Treatment) %>%
  summarize(mean = mean(Corrected_TA),
            sd = sd(Corrected_TA, ))


TANight<-Night %>%
  group_by(Month, Treatment) %>%
  summarize(mean = mean(Corrected_TA),
            sd = sd(Corrected_TA))

####################--------------------------------------------------------------
#Salinity
####################--------------------------------------------------------------
Day %>%
  group_by(Month, Treatment) %>%
  summarize(mean = mean(Salinity),
            sd = sd(Salinity))

Night %>%
  group_by(Month, Treatment) %>%
  summarize(mean = mean(Salinity),
            sd = sd(Salinity))

####################--------------------------------------------------------------
#Dissolved Oxygen
####################--------------------------------------------------------------
Day %>%
  group_by(Month, Treatment) %>%
  summarize(mean = mean(DO),
            sd = sd(DO))

Night %>%
  group_by(Month, Treatment) %>%
  summarize(mean = mean(DO),
            sd = sd(DO))


####################--------------------------------------------------------------
#Carbonate Chemistry | DIC, pCO2, Calcite, and Aragonite | Projected Values
####################--------------------------------------------------------------
pH_Meta<- read.csv("pH_Meta.csv")

Day<-filter(pH_Meta, DayNight == "Day" )
Night<-filter(pH_Meta, DayNight == "Night")


####################--------------------------------------------------------------
#TIMEPOINT 3 for each Month and Treatment
####################--------------------------------------------------------------

####################--------------------------------------------------------------
#Day Pool DIC + Error
####################--------------------------------------------------------------

DayPoolCO2<-carb(flag=8, Day$HannapH_Tot, Day$Corrected_TA/1000000, S=Day$Salinity, T=Day$Temp, Patm=1, P=0, 
                 Pt=Day$Phosphate/1000000, Sit=0, k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")
#TA is divided by 1000 because all calculations are in mol/kg in the seacarb package

# calculate error propagation
Dayer<-errors(flag=8, Day$HannapH_Tot, Day$Corrected_TA/1000000, 
              S=Day$Salinity,  T=Day$Temp, Patm=1, P=0, 
              Pt=Day$Phosphate/1000000, Sit=0, evar1 = 0.01, evar2 = 5e-6) 

#average error for DIC based on pH and TA
mean(Dayer$DIC*1000000)

sd(Dayer$DIC*1000000)/sqrt(nrow(Dayer))


#convert CO2, HCO3, CO3, DIC, and Alk back to micromol for easier interpretation
DayPoolCO2[,c("CO2","HCO3","CO3","DIC","ALK")]<-DayPoolCO2[,c("CO2","HCO3","CO3","DIC","ALK")]*1000000

#combine all the data into new data sheet
DayCarbChem<-cbind(Day,DayPoolCO2[-c(1:6,11:13)])

write.csv(DayCarbChem, file="DayCarbChem_Mussel.csv")

####################--------------------------------------------------------------
#Night Pool DIC + Error
####################--------------------------------------------------------------

#Night Pool DIC
NightPoolCO2<-carb(flag=8, Night$HannapH, Night$CorrectedTA/1000000, S=Night$Salinity, T=Night$Temp, Patm=1, P=0, 
                   Pt=Night$Phosphate/1000000, Sit=0, k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")

# calculate error propagation
Nighter<-errors(flag=8, Night$HannapH, Night$CorrectedTA/1000000, 
                S=Night$Salinity,  T=Night$Temp, Patm=1, P=0, 
                Pt=Night$Phosphate/1000000, Sit=0, evar1 = 0.01, evar2 = 5e-6) 


#average error for DIC based on pH and TA
mean(Nighter$DIC*1000000)

sd(Nighter$DIC*1000000)/sqrt(nrow(Nighter))


#convert CO2, HCO3, CO3, DIC, and Alk back to micromol for easier interpretation
NightPoolCO2[,c("CO2","HCO3","CO3","DIC","ALK")]<-NightPoolCO2[,c("CO2","HCO3","CO3","DIC","ALK")]*1000000

NightCarbChem<-cbind(Night,NightPoolCO2[-c(1:6,11:13)])

#combine all the data into new data sheet
write.csv(NightCarbChem, file="NightCarbChem_Mussel.csv")

####################--------------------------------------------------------------
#TIMEPOINT 1 - For Future "Projected pH values" Mussels | From Timepoint 1
####################--------------------------------------------------------------
WaterMeta<- read.csv("WaterMeta_DIC_Musselsl.csv")

Day<-filter(WaterMeta, DayNight == "Day" )
Night<-filter(WaterMeta, DayNight == "Night")

####################--------------------------------------------------------------
#Day Pool DIC + Error Timepoint 1
####################--------------------------------------------------------------
#Calculate DIC first
DayPoolCO2<-carb(flag=8, Day$HannapH, Day$CorrectedTA/1000000, S=Day$Salinity, T=Day$Temp, Patm=1, P=0, 
                 Pt=Day$Phosphate/1000000, Sit=0, k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")
#TA is divided by 1000 because all calculations are in mol/kg in the seacarb package

# calculate error propagation
Dayer<-errors(flag=8, Day$HannapH, Day$CorrectedTA/1000000, 
              S=Day$Salinity,  T=Day$Temp, Patm=1, P=0, 
              Pt=Day$Phosphate/1000000, Sit=0, evar1 = 0.01, evar2 = 5e-6) 

#average error for DIC based on pH and TA
mean(Dayer$DIC*1000000)
#8.383398
sd(Dayer$DIC*1000000)/sqrt(nrow(Dayer))
#0.8598083

#convert CO2, HCO3, CO3, DIC, and Alk back to micromol for easier interpretation
DayPoolCO2[,c("CO2","HCO3","CO3","DIC","ALK")]<-DayPoolCO2[,c("CO2","HCO3","CO3","DIC","ALK")]*1000000

#average DIC based on pH and TA for Timepoint 1
mean(DayPoolCO2$DIC)
#1848.601

#combine all the data
DayCarbChem<-cbind(Day,DayPoolCO2[-c(1:6,11:13)])

####################--------------------------------------------------------------
#Day Pool DIC PROJECTED
####################--------------------------------------------------------------

DayPoolCO2<-carb(flag=8, Day$HannapH, Day$CorrectedTA/1000000, S=Day$Salinity, T=Day$Temp, Patm=1, P=0, 
                 Pt=Day$Phosphate/1000000, Sit=0, k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")
#TA is divided by 1000000 because all calculations are in mol/kg in the seacarb package

#average error for DIC based on pH and TA
#mean(DayPoolCO2$DIC*1000000)

#8.383398
sd(Dayer$DIC*1000000)/sqrt(nrow(Dayer))
#0.8598083

#Add 200 micromol so 0.0002 mol/kg to DIC
DayPoolCO2$DICTwo<-DayPoolCO2[, 16] + 0.0002

#Back Calculate pH from DIC and pCO2 - DIC in mol/kg
DayPoolCO2.pH<-carb(flag=15, DayPoolCO2$DICTwo, Day$CorrectedTA/1000000, S=Day$Salinity, T=Day$Temp, Patm=1, P=0, 
                    Pt=Day$Phosphate/1000000, Sit=0, k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")
#TA is divided by 1000 because all calculations are in mol/kg in the seacarb package

#convert CO2, HCO3, CO3, DIC, and Alk back to micromol for easier interpretation
DayPoolCO2.pH[,c("CO2","HCO3","CO3","DIC","ALK")]<-DayPoolCO2.pH[,c("CO2","HCO3","CO3","DIC","ALK")]*1000000

#average error for DIC based on pH and TA
mean(DayPoolCO2.pH$pH)
#7.616471

#combine all the data
DayCarbChem.newpH<-cbind(Day,DayPoolCO2.pH[-c(1:5,9:13)])

write.csv(DayCarbChem.newpH, file="DayCarbChemProjected.csv")


####################--------------------------------------------------------------
#Night Pool DIC + Error Timepoint 1
####################--------------------------------------------------------------
#Calculate DIC first
NightPoolCO2<-carb(flag=8, Night$HannapH, Night$CorrectedTA/1000000, S=Night$Salinity, T=Night$Temp, Patm=1, P=0, 
                   Pt=Night$Phosphate/1000000, Sit=0, k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")
#TA is divided by 1000 because all calculations are in mol/kg in the seacarb package

# calculate error propagation
Nighter<-errors(flag=8, Night$HannapH, Night$CorrectedTA/1000000, 
                S=Night$Salinity,  T=Night$Temp, Patm=1, P=0, 
                Pt=Night$Phosphate/1000000, Sit=0, evar1 = 0.01, evar2 = 5e-6) 

#average error for DIC based on pH and TA
mean(Nighter$DIC*1000000)


#convert CO2, HCO3, CO3, DIC, and Alk back to micromol for easier interpretation
NightPoolCO2[,c("CO2","HCO3","CO3","DIC","ALK")]<-NightPoolCO2[,c("CO2","HCO3","CO3","DIC","ALK")]*1000000

#average DIC based on pH and TA for Timepoint 1
mean(NightPoolCO2$DIC)
#1968.943

#combine all the data
NightCarbChem<-cbind(Night,NightPoolCO2[-c(1:6,11:13)])


####################--------------------------------------------------------------
#Night Pool DIC PROJECTED
####################--------------------------------------------------------------

NightPoolCO2<-carb(flag=8, Night$HannapH, Night$CorrectedTA/1000000, S=Night$Salinity, T=Night$Temp, Patm=1, P=0, 
                   Pt=Night$Phosphate/1000000, Sit=0, k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")
#TA is divided by 1000000 because all calculations are in mol/kg in the seacarb package

#average error for DIC based on pH and TA
mean(NightPoolCO2$DIC*1000000)
#1968.943


#Add 200 micromol so 0.0002 mol/kg to DIC
NightPoolCO2$DICTwo<-NightPoolCO2[, 16] + 0.0002

#Back Calculate pH from DIC and pCO2 - DIC in mol/kg
NightPoolCO2.pH<-carb(flag=15, NightPoolCO2$DICTwo, Night$CorrectedTA/1000000, S=Night$Salinity, T=Night$Temp, Patm=1, P=0, 
                      Pt=Night$Phosphate/1000000, Sit=0, k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")
#TA is divided by 1000 because all calculations are in mol/kg in the seacarb package

#convert CO2, HCO3, CO3, DIC, and Alk back to micromol for easier interpretation
NightPoolCO2.pH[,c("CO2","HCO3","CO3","DIC","ALK")]<-NightPoolCO2.pH[,c("CO2","HCO3","CO3","DIC","ALK")]*1000000

#average error for DIC based on pH and TA
mean(NightPoolCO2.pH$pH)
#7.866437

#combine all the data
NightCarbChem.newpH<-cbind(Night,NightPoolCO2.pH[-c(1:5,9:13)])

write.csv(NightCarbChem.newpH, file="NightCarbChemProjected.csv")


####################--------------------------------------------------------------
#Day Pool DIC PROJECTED - Text Values
####################--------------------------------------------------------------
DayProject<- read.csv("DayCarbChemProjected.csv")

pH.36<-filter(DayProject, !Month =="March") # Only comparing when manipulations were on.


#CO2 pool difference between no CO2 - Daytime
pH.36 %>%
  group_by(CO2) %>%
  summarize(mean = mean(pH),
            sd = sd(pH),
            n = n(),
            se=sd/sqrt(n))

pH.36 %>%
  group_by(CO2) %>%
  summarize(mean = mean(DIC),
            sd = sd(DIC),
            n = n(),
            se=sd/sqrt(n),
            DIC_range_min = min(DIC),
            DIC_range_max = max(DIC))

####################--------------------------------------------------------------
#CO2 pool difference between no CO2 - Nighttime
####################--------------------------------------------------------------

NightProject<- read.csv("NightCarbChemProjected.csv")

pH.36Night<-filter(NightProject, !Month =="March")

pH.36Night %>%
  group_by(CO2) %>%
  summarize(mean = mean(pH),
            sd = sd(pH),
            n = n(),
            se=sd/sqrt(n))

pH.36Night %>%
  group_by(CO2) %>%
  summarize(mean = mean(DIC),
            sd = sd(DIC),
            n = n(),
            se=sd/sqrt(n),
            DIC_range_min = min(DIC),
            DIC_range_max = max(DIC))

####################--------------------------------------------------------------
#DayTime
####################--------------------------------------------------------------
DayCarbChem<- read.csv("DayCarbChem_Mussel.csv")

DayCarbChem$Treatment<-factor(DayCarbChem$Treatment, levels = c("Unmanipulated", "CO2 Addition","Warming", "Warming + CO2"), labels=c("Unmanipulated", "CO2 Added","Warmed", "CO2 + Warmed"))

DayCarbChem$Month<-factor(DayCarbChem$Month, levels = c("March","May", "June", "July", "August", "September"), labels = c("March","May", "June", "July", "Aug.", "Sept."))

DayDIC<- DayCarbChem %>%
  group_by(Month, Treatment) %>%
  summarize(mean = round(mean(DIC), 0),
            sd = round(sd(DIC), 0),
            n=n())


DaypCO2<- DayCarbChem %>%
  group_by(Month, Treatment) %>%
  summarize(mean = round(mean(pCO2), 2),
            sd = round(sd(pCO2), 2),
            n=n())

DayArg<- DayCarbChem %>%
  group_by(Month, Treatment) %>%
  summarize(mean = round(mean(OmegaAragonite), 2),
            sd = round(sd(OmegaAragonite), 2),
            n=n())

DayCal<- DayCarbChem %>%
  group_by(Month, Treatment) %>%
  summarize(mean = round(mean(OmegaCalcite), 2),
            sd = round(sd(OmegaCalcite), 2),
            n=n())

####################--------------------------------------------------------------
#Nighttime
####################--------------------------------------------------------------
NightCarbChem<- read.csv("NightCarbChem_Mussel.csv")

NightCarbChem$Treatment<-factor(NightCarbChem$Treatment, levels = c("Unmanipulated", "CO2 Addition","Warming", "Warming + CO2"), labels=c("Unmanipulated", "CO2 Added","Warmed", "CO2 + Warmed"))

NightCarbChem$Month<-factor(NightCarbChem$Month, levels = c("March", "July", "September"))

NightDIC<- NightCarbChem %>%
  group_by(Month, Treatment) %>%
  summarize(mean = round(mean(DIC), 0),
            sd = round(sd(DIC), 0),
            n=n())

NightpCO2<- NightCarbChem %>%
  group_by(Month, Treatment) %>%
  summarize(mean = round(mean(pCO2), 0),
            sd = round(sd(pCO2), 0),
            n=n())


NightArg<- NightCarbChem %>%
  group_by(Month, Treatment) %>%
  summarize(mean = round(mean(OmegaAragonite), 2),
            sd = round(sd(OmegaAragonite), 2),
            n=n())

NightCal<- NightCarbChem %>%
  group_by(Month, Treatment) %>%
  summarize(mean = round(mean(OmegaCalcite), 2),
            sd = round(sd(OmegaCalcite), 2),
            n=n())


####################--------------------------------------------------------------
#pH Data Analysis
####################--------------------------------------------------------------
pH.Zero<-filter(pH_Meta, Month == "March")
pH.Three<-filter(pH_Meta, Month =="July")
pH.Six<-filter(pH_Meta, Month =="September")

#Check distribution
hist(pH.Zero$HannapH_Tot)
shapiro.test(pH.Zero$HannapH_Tot) # Normal p-value = 0.06727

hist(pH.Three$HannapH_Tot)
shapiro.test(pH.Three$HannapH_Tot) # Not normal p-value = 0.005146 #GAMMA

hist(pH.Six$HannapH_Tot)
shapiro.test(pH.Six$HannapH_Tot) # Normal p-value = 0.08593

####################--------------------------------------------------------------
#GLM with Time 0 
####################--------------------------------------------------------------
M.ph.Zero<-glm(HannapH_Tot~CO2+Heat+DayNight+Heat*CO2+Heat*DayNight+CO2*DayNight+Volume, data=pH.Zero)
summary(M.ph.Zero)
anova(M.ph.Zero, test="F")
#DayNight Sig - p<0.001

check_model(M.ph.Zero) #Residuals fit

M.phZero.post<- emmeans(M.ph.Zero, ~CO2*Heat|DayNight)

m_tukey=contrast(M.phZero.post, method='pairwise', by=NULL)
summary(m_tukey)

simulationOutput<-simulateResiduals(fittedModel = M.ph.Three)
plot(simulationOutput)


####################--------------------------------------------------------------
#GLM with Time 3 
####################--------------------------------------------------------------
M.ph.Three<-glm(HannapH_Tot~CO2+Heat+DayNight+Heat*CO2+Heat*DayNight+CO2*DayNight+Volume, family=Gamma(link='log'), data=pH.Three)
summary(M.ph.Three)
anova(M.ph.Three,test="F")
#DayNight + CO2 p<0.001, Heat p=0.090, Heat*CO2 p=0.0428, CO2*DN p=0.009

check_model(M.ph.Three) #Residuals fit 

M.phThree.post<- emmeans(M.ph.Three, ~CO2*Heat|DayNight)

m_tukey=contrast(M.phThree.post, method='pairwise', by=NULL)
summary(m_tukey)


####################--------------------------------------------------------------
#GLM with Time 6 
####################--------------------------------------------------------------
M.ph.Six<-glm(HannapH_Tot~CO2+Heat+DayNight+Heat*CO2+Heat*DayNight+CO2*DayNight+Volume, data=pH.Six)
summary(M.ph.Six)
anova(M.ph.Six, test="F")
#DayNight p<0.001
check_model(M.ph.Six) #Residuals fit 

M.phSix.post<- emmeans(M.ph.Six, ~CO2*Heat|DayNight) 

m_tukey=contrast(M.phSix.post, method='pairwise', by=NULL)
summary(m_tukey)


####################--------------------------------------------------------------
#TEMPERATURE DATA
####################--------------------------------------------------------------

####################--------------------------------------------------------------
#Temperature Data Analysis
####################--------------------------------------------------------------
Temp<-read.csv("DetidedTemp.csv")
####################--------------------------------------------------------------
#Temperature Data Analysis
####################--------------------------------------------------------------
Day<-filter(Temp, DayNight == "Day" )
Night<-filter(Temp, DayNight == "Night")


T.Zero<-filter(Temp, Timepoint == "0")
T.Three<-filter(Temp, Timepoint =="3")
T.Six<-filter(Temp, Timepoint =="6")

#Without March for Daytime Tides after manipulations are on
T.36D<-filter(Day, !Timepoint=="0")
T.36N<-filter(Night, !Timepoint =="0")

#Without March for All Tides after manipulations are on
T.36All<-filter(Temp, !Timepoint=="0")

#Looking at just Time 0 / March - Pre-manipulation
T.Zero %>%
  group_by(Heat) %>%
  summarize(mean = mean(Daily99MaxPercentile),
            sd = sd(Daily99MaxPercentile))

#Calculating difference between Time 3 (July) and Time 6 (September) after manipulations are on - only for Daytime Tides
T.36D %>%
  group_by(Heat) %>%
  summarize(mean = mean(Daily99MaxPercentile),
            sd = sd(Daily99MaxPercentile),
            n = n(),
            se=sd/sqrt(n))

t.test(Daily99MaxPercentile~Heat, data=T.36D) # No sig. differences

#Calculating difference between Time 3 (July) and Time 6 (September) after manipulations are on - only for Nighttime Tides
T.36N %>%
  group_by(Heat) %>%
  summarize(mean = mean(Daily99MaxPercentile),
            sd = sd(Daily99MaxPercentile),
            n = n(),
            se=sd/sqrt(n))

#Calculating difference between Time 3 (July) and Time 6 (September) after manipulations are on - for All Tides
T.36All %>%
  group_by(Heat) %>%
  summarize(mean = mean(Daily99MaxPercentile),
            sd = sd(Daily99MaxPercentile),
            n = n(),
            se=sd/sqrt(n))
#All day tides including timepoint 0 (March)
Day %>%
  group_by(Heat) %>%
  summarize(mean = mean(Daily99MaxPercentile),
            sd = sd(Daily99MaxPercentile),
            n = n(),
            se=sd/sqrt(n))

#All nighttime tides including timepoint 0 (March)
Night %>%
  group_by(Heat) %>%
  summarize(mean = mean(Daily99MaxPercentile),
            sd = sd(Daily99MaxPercentile),
            n = n(),
            se=sd/sqrt(n))

#Check distribution
hist(T.Zero$Daily99MaxPercentile)
shapiro.test(T.Zero$Daily99MaxPercentile) # not normal  p-value = 0.02642, GAMMA

hist(T.Three$Daily99MaxPercentile)
shapiro.test(T.Three$Daily99MaxPercentile) #normal p-value =0.05309

hist(T.Six$Daily99MaxPercentile)
shapiro.test(T.Six$Daily99MaxPercentile) #not normal p-value = 0.005003

####################--------------------------------------------------------------
#GLM with Time 0 
####################--------------------------------------------------------------
M.T.Zero<-glm(Daily99MaxPercentile~CO2+Heat+DayNight+Heat*CO2+Heat*DayNight+CO2*DayNight+Volume, family=Gamma(link='log'), data=T.Zero)
summary(M.T.Zero)
anova(M.T.Zero, test="F")
#No difference beside DayNight p<0.001

check_model(M.T.Zero) #Pretty okay...

M1.T.Zero.post<- emmeans(M.T.Zero, ~CO2*Heat|DayNight)

m_tukey=contrast(M1.T.Zero.post, method='pairwise', by=NULL)
summary(m_tukey)


####################--------------------------------------------------------------
#GLM with Time 3 
####################--------------------------------------------------------------
M.T.Three<-glm(Daily99MaxPercentile~CO2+Heat+DayNight+Heat*CO2+Heat*DayNight+CO2*DayNight+Volume,family=Gamma(link='log'), data=T.Three)
summary(M.T.Three)
anova(M.T.Three, test="F")
#No difference beside DayNight p<0.001

check_model(M.T.Three) #Decent

M1.T.Three.post<- emmeans(M.T.Three, ~CO2*Heat|DayNight)

m_tukey=contrast(M1.T.Three.post, method='pairwise', by=NULL)
summary(m_tukey)

####################--------------------------------------------------------------
#GLM with Time 6 
####################--------------------------------------------------------------
M.T.Six<-glm(Daily99MaxPercentile~CO2+Heat+DayNight+Heat*CO2+Heat*DayNight+CO2*DayNight+Volume,family=Gamma(link='log'), data=T.Six)

summary(M.T.Six)
anova(M.T.Six, test="F")


check_model(M.T.Six2) #Decent

M1.T.Six.post<- emmeans(M.T.Six, ~CO2*Heat|DayNight)

m_tukey=contrast(M1.T.Six.post, method='pairwise', by=NULL)
summary(m_tukey)


####################--------------------------------------------------------------
#Temp without any other factors
####################--------------------------------------------------------------
#Check distribution
hist(T.Zero$Daily99MaxPercentile)
shapiro.test(T.Zero$Daily99MaxPercentile) #not normal  p-value = 0.02642 GAMMA

hist(T.Three$Daily99MaxPercentile)
shapiro.test(T.Three$Daily99MaxPercentile) # Close to normal p-value = 0.05309 

hist(T.Six$Daily99MaxPercentile)
shapiro.test(T.Six$Daily99MaxPercentile) #not normal p-value = 0.005003 GAMMA

####################--------------------------------------------------------------
#GLM with Time 0 
####################--------------------------------------------------------------
M.T.Zero<-glm(Daily99MaxPercentile~Heat, family=Gamma(link='log'), data=T.Zero)
summary(M.T.Zero)
anova(M.T.Zero, test="F")
#No difference

check_model(M.T.Zero) #Pretty okay...



####################--------------------------------------------------------------
#GLM with Time 3 
####################--------------------------------------------------------------
M.T.Three<-glm(Daily99MaxPercentile~Heat, data=T.Three)
summary(M.T.Three)
anova(M.T.Three, test="F")
#No difference

check_model(M.T.Three) #Decent


####################--------------------------------------------------------------
#GLM with Time 6 
####################--------------------------------------------------------------
M.T.Six<-glm(Daily99MaxPercentile~Heat,family=Gamma(link='log'), data=T.Six)

summary(M.T.Six)
anova(M.T.Six, test="F")
#No difference

check_model(M.T.Six) #Decent

M1.T.Six.post<- emmeans(M.T.Six, ~CO2*Heat|DayNight)

m_tukey=contrast(M1.T.Six.post, method='pairwise', by=NULL)
summary(m_tukey)




####################--------------------------------------------------------------
#Average Temperature Data Analysis
####################--------------------------------------------------------------
#Check distribution
hist(T.Zero$DailyAverage)
shapiro.test(T.Zero$DailyAverage) #normal  p-value = 0.1267

hist(T.Three$DailyAverage)
shapiro.test(T.Three$DailyAverage) # normal p-value = 0.06715 

hist(T.Six$DailyAverage)
shapiro.test(T.Six$DailyAverage) #not normal p-value = 0.02063 GAMMA


####################--------------------------------------------------------------
#GLM with Time 0 
####################--------------------------------------------------------------
M.T.Zero<-glm(DailyAverage~CO2+Heat+DayNight+Heat*CO2+Heat*DayNight+CO2*DayNight, data=T.Zero)
summary(M.T.Zero)
anova(M.T.Zero, test="F")
#No difference beside DayNight p<0.001

check_model(M.T.Zero) #Pretty okay...

M1.T.Zero.post<- emmeans(M.T.Zero, ~CO2*Heat|DayNight)

m_tukey=contrast(M1.T.Zero.post, method='pairwise', by=NULL)
summary(m_tukey)


####################--------------------------------------------------------------
#GLM with Time 3 
####################--------------------------------------------------------------
M.T.Three<-glm(DailyAverage~CO2+Heat+DayNight+Heat*CO2+Heat*DayNight+CO2*DayNight, data=T.Three)
summary(M.T.Three)
anova(M.T.Three, test="F")
#No difference beside DayNight p<0.001

check_model(M.T.Three) #Decent

M1.T.Three.post<- emmeans(M.T.Three, ~CO2*Heat|DayNight)

m_tukey=contrast(M1.T.Three.post, method='pairwise', by=NULL)
summary(m_tukey)

####################--------------------------------------------------------------
#GLM with Time 6 
####################--------------------------------------------------------------
M.T.Six<-glm(DailyAverage~CO2+Heat+DayNight+Heat*CO2+Heat*DayNight+CO2*DayNight,family=Gamma(link='log'), data=T.Six)

summary(M.T.Six)
anova(M.T.Six, test="F") # Heat - p=0.03411

check_model(M.T.Six) #Decent

M1.T.Six.post<- emmeans(M.T.Six, ~CO2*Heat|DayNight)

m_tukey=contrast(M1.T.Six.post, method='pairwise', by=NULL)
summary(m_tukey)


####################--------------------------------------------------------------
#Shell Strength Data
####################--------------------------------------------------------------
setwd("/Users/racinerangel/Desktop/Climate Manipulations/Shell analyses")
####################--------------------------------------------------------------
Muss_Str <- read.csv("ShellStrength_CL.csv")

#DO NOT FORGET***************Reorder treatments
Muss_Str$Treatment<-factor(Muss_Str$Treatment, levels = c("Control", "Heat", "CO2", "Both"), labels = c("Unmanipulated", "Warming", "CO2 Addition", "Warming + CO2"))

####################--------------------------------------------------------------
#Average Strength
####################--------------------------------------------------------------
Strength.avg <- Muss_Str %>%
  group_by(Month,Timepoint, Treatment) %>%
  summarize(mean = mean(Strength.Stand),
            n = sum(Strength.Stand),
            sd= sd(Strength.Stand),
            se=sd/sqrt(n))


####################--------------------------------------------------------------
#GLM with Time 0 and then 3 vs 6 Time of exposure
####################--------------------------------------------------------------
Str.Zero<-filter(Muss_Str, Timepoint == "0")
Str.Three<-filter(Muss_Str, Timepoint =="3")
Str.Six<-filter(Muss_Str, Timepoint =="6")

#Check distribution
hist(Str.Zero$Strength.Stand)
shapiro.test(Str.Zero$Strength.Stand) #normal  p-value = 0.1549

hist(Str.Three$Strength.Stand)
shapiro.test(Str.Three$Strength.Stand) #normal p-value = 0.1989

hist(Str.Six$Strength.Stand)
shapiro.test(Str.Six$Strength.Stand) #normal  p-value = 0.8741

####################--------------------------------------------------------------
#GLM with Time 0 
####################--------------------------------------------------------------
M.Str.Zero<-glm(Strength.Stand~CO2*Heat+Volume, data=Str.Zero)
summary(M.Str.Zero)
anova(M.Str.Zero, test="F")
#No difference

quartz()
check_model(M.Str.Zero) #

M1.Str.Zero.post<- emmeans(M.Str.Zero, ~CO2*Heat)

m_tukey=contrast(M1.Str.Zero.post, method='pairwise', by=NULL)
summary(m_tukey)

####################--------------------------------------------------------------
#GLM with Time 3 
####################--------------------------------------------------------------
M.Str.Three<-glm(Strength.Stand~CO2+Heat+Heat*CO2+Volume, data=Str.Three)
summary(M.Str.Three)
anova(M.Str.Three, test="F")
#No difference

quartz()
check_model(M.Str.Three) #variance a little off

M1.Str.Three.post<- emmeans(M.Str.Three, ~CO2*Heat)

m_tukey=contrast(M1.Str.Three.post, method='pairwise', by=NULL)
summary(m_tukey)

####################--------------------------------------------------------------
#GLM with Time 6 
####################--------------------------------------------------------------
M.Str.Six <- glm(Strength.Stand~CO2*Heat+Volume, data=Str.Six)
summary(M.Str.Six)
anova(M.Str.Six, test="F")
#Heat sig: p=0.029757, CO2*Heat: p=0.004272


quartz()
check_model(M.Str.Six) #decent but not great

M1.Str.Six.post<- emmeans(M.Str.Six, ~CO2*Heat)

m_tukey=contrast(M1.Str.Six.post, method='pairwise', by=NULL)
summary(m_tukey)

####################--------------------------------------------------------------
#PERCENT DIFFERENCE BETWEEN GROUPS - need library(scales)
####################--------------------------------------------------------------
Str.Six %>%
  group_by(Treatment) %>%
  summarise(mean_str = mean(Strength.Stand)) %>%
  mutate(perc_change = scales::percent((mean_str[Treatment == "CO2 Addition"] - mean_str[Treatment == "Warming + CO2"])/mean_str[Treatment == "CO2 Addition"]))

####################--------------------------------------------------------------
#SHELL THICKNESS
####################--------------------------------------------------------------
setwd("/Users/racinerangel/Desktop/Climate Manipulations/Shell Analyses/Shell Thickness")
####################--------------------------------------------------------------

ShellMeta<-read.csv("Sitka_Mussel_CLShells.csv")
ShellMeta$Timepoint<-as.factor(ShellMeta$Timepoint)

#ShellMeta$Month<-factor(ShellMeta$Month, levels = c("March", "July", "Sept"), labels = c("0", "3", "6")) #For plotting

ShellMeta$Treatment<-factor(ShellMeta$Treatment, levels = c("Control", "Heat", "CO2", "Both"), labels=c("Unmanipulated", "Warming", "CO2 Addition", "Warming + CO2"))

####################--------------------------------------------------------------
#Shell thickness Summary Stats
####################--------------------------------------------------------------

ShellTH_L <- ShellMeta %>%
  group_by(Month, Timepoint, Treatment) %>%
  summarize(mean = mean(StandardA),
            sd = sd(StandardA),
            n = n(),
            se = sd/sqrt(n))

ShellTH_M <- ShellMeta %>%
  group_by(Month, Timepoint, Treatment) %>%
  summarize(mean = mean(StandardB),
            sd = sd(StandardB),
            n = n(),
            se = sd/sqrt(n))

ShellTH_C <- ShellMeta %>%
  group_by(Month, Timepoint, Treatment) %>%
  summarize(mean = mean(StandardC),
            sd = sd(StandardC),
            n = n(),
            se = sd/sqrt(n))

ShellTH_TL <- ShellMeta %>%
  group_by(Month, Timepoint, Treatment) %>%
  summarize(mean = round(mean(Total_L), 3),
            sd = round(sd(Total_L), 3),
            n = n(),
            se = round(sd/sqrt(n), 3))

ShellTH_TLTreat <- ShellMeta %>%
  group_by(Treatment) %>%
  summarize(mean = round(mean(Total_L), 3),
            sd = round(sd(Total_L), 3),
            n = n(),
            se = round(sd/sqrt(n), 3))

####################--------------------------------------------------------------
#GLM with Time 0 and then 3 vs 6 Time of exposure
####################--------------------------------------------------------------
Lip.Zero<-filter(ShellMeta, Month == "March")
Lip.Three<-filter(ShellMeta, Month =="July")
Lip.Six<-filter(ShellMeta, Month =="September")

#Check distribution
hist(Lip.Zero$StandardA)
shapiro.test(Lip.Zero$StandardA) #normal p-value = 0.7997

hist(Lip.Three$StandardA)
shapiro.test(Lip.Three$StandardA) #normal p-value = 0.1323

hist(Lip.Six$StandardA)
shapiro.test(Lip.Six$StandardA) #not normal p-value = 0.002466 GAMMA

#Lip.Six$Tr.StandardA<-log10(Lip.Six$StandardA) #Transformed 
#shapiro.test(Lip.Six$Tr.StandardA) #normal now
#hist(Lip.Six$Tr.StandardA)

####################--------------------------------------------------------------
#GLM with Time 0 
####################--------------------------------------------------------------
M.Lip.Zero<-glm(StandardA~CO2*Heat+Volume, data=Lip.Zero)
summary(M.Lip.Zero)
anova(M.Lip.Zero, test="F")
# no difference


check_model(M.Lip.Zero) #Goodish

M1.Lip.Zero.post<- emmeans(M.Lip.Zero, ~CO2*Heat)

m_tukey=contrast(M1.Lip.Zero.post, method='pairwise', by=NULL)
summary(m_tukey)


####################--------------------------------------------------------------
#GLM with Time 3 
####################--------------------------------------------------------------
M.Lip.Three<-glm(StandardA~CO2*Heat+Volume, data=Lip.Three)
summary(M.Lip.Three)
anova(M.Lip.Three, test="F")


check_model(M.Lip.Three)

#CO2 * Heat Sig Interaction

M1.Lip.Three.post<- emmeans(M.Lip.Three, ~CO2*Heat)

m_tukey=contrast(M1.Lip.Three.post, method='pairwise', by=NULL)
summary(m_tukey)

####################--------------------------------------------------------------
#GLM with Time 6 
####################--------------------------------------------------------------
M.Lip.Six <- glm(StandardA~CO2*Heat+Volume, family=Gamma(link='log'), data=Lip.Six)
summary(M.Lip.Six)
anova(M.Lip.Six, test="F")


check_model(M.Lip.Six)

M1.Lip.Six.post<- emmeans(M.Lip.Six, ~CO2*Heat)

m_tukey=contrast(M1.Lip.Six.post, method='pairwise', by=NULL)
summary(m_tukey)

####################--------------------------------------------------------------
#MIDDLE SHELL THICKNESS
####################--------------------------------------------------------------
Mid.Zero<-filter(ShellMeta, Month == "March")
Mid.Three<-filter(ShellMeta, Month == "July")
Mid.Six<-filter(ShellMeta, Month == "September")

#Check distribution
hist(Mid.Zero$StandardB)
shapiro.test(Mid.Zero$StandardB) #normal p-value = 0.12

hist(Mid.Three$StandardB)
shapiro.test(Mid.Three$StandardB) #normal p-value = 0.4753

hist(Mid.Six$StandardB)
shapiro.test(Mid.Six$StandardB) #normal p-value = 0.9625

####################--------------------------------------------------------------
#GLM with Time 0 
####################--------------------------------------------------------------
M.Mid.Zero<-glm(StandardB~CO2+Heat+Heat*CO2+Volume, data=Mid.Zero)
summary(M.Mid.Zero)
anova(M.Mid.Zero, test="F")
# no difference


check_model(M.Mid.Zero)

M1.Mid.Zero.post<- emmeans(M.Mid.Zero, ~CO2*Heat)

m_tukey=contrast(M1.Mid.Zero.post, method='pairwise', by=NULL)
summary(m_tukey) #CO2 vs Unmanpulated Sig Different

####################--------------------------------------------------------------
#GLM with Time 3 
####################--------------------------------------------------------------
M.Mid.Three<-glm(StandardB~CO2*Heat+Volume, data=Mid.Three)
summary(M.Mid.Three)
anova(M.Mid.Three, test="F")
#CO2 Sig

check_model(M.Mid.Three)

fm1 <- aov(StandardB~Treatment, data=Mid.Three)
anova(fm1)

M1.Mid.Three.post<- emmeans(M.Mid.Three, ~CO2*Heat)

m_tukey=contrast(M1.Mid.Three.post, method='pairwise', by=NULL)

summary(m_tukey) #CO2 vs Unmanpulated almost Sig Different

####################--------------------------------------------------------------
#GLM with Time 6 
####################--------------------------------------------------------------
M.Mid.Six <- glm(StandardB~CO2*Heat+Volume, data=Mid.Six)
summary(M.Mid.Six)
anova(M.Mid.Six, test="F")

fm1 <- aov(StandardB~Treatment, data=Mid.Six)
anova(fm1)

#CO2 Sig

check_model(M.Mid.Six)

M1.Mid.Six.post<- emmeans(M.Mid.Six, ~CO2*Heat)

m_tukey=contrast(M1.Mid.Six.post, method='pairwise', by=NULL)
summary(m_tukey) #CO2 vs Unmanpulated Sig Different

####################--------------------------------------------------------------
#BASE SHELL THICKNESS
####################--------------------------------------------------------------
Base.Zero<-filter(ShellMeta, Month == "March")
Base.Three<-filter(ShellMeta, Month =="July")
Base.Six<-filter(ShellMeta, Month =="September")

#Check distribution
hist(Base.Zero$StandardC)
shapiro.test(Base.Zero$StandardC) #normal p-value = 0.4872

hist(Base.Three$StandardC)
shapiro.test(Base.Three$StandardC) #normal p-value = 0.1852

hist(Base.Six$StandardC)
shapiro.test(Base.Six$StandardC) # normal p-value = 0.6154

####################--------------------------------------------------------------
#GLM with Time 0 
####################--------------------------------------------------------------
M.Base.Zero<-glm(StandardC~CO2*Heat+Volume, data=Base.Zero)
summary(M.Base.Zero)
anova(M.Base.Zero, test="F")
#No difference

check_model(M.Base.Zero)

M1.Base.Zero.post<- emmeans(M.Base.Zero, ~CO2*Heat)

m_tukey=contrast(M1.Base.Zero.post, method='pairwise', by=NULL)
summary(m_tukey) 

####################--------------------------------------------------------------
#GLM with Time 3 
####################--------------------------------------------------------------
M.Base.Three<-glm(StandardC~CO2*Heat+Volume, data=Base.Three)
summary(M.Base.Three)
anova(M.Base.Three, test="F")
#No difference

check_model(M.Base.Three) #homogeneity a little curved

M1.Base.Three.post<- emmeans(M.Base.Three, ~CO2*Heat)

m_tukey=contrast(M1.Base.Three.post, method='pairwise', by=NULL)
summary(m_tukey) 

####################--------------------------------------------------------------
#GLM with Time 6 
####################--------------------------------------------------------------
M.Base.Six <- glm(StandardC~CO2*Heat+Volume, data=Base.Six)
summary(M.Base.Six)
anova(M.Base.Six, test="F")
#CO2 and Heat marg Sig

check_model(M.Base.Six) #Decreases a bit (homogeneity)

M1.Base.Six.post<- emmeans(M.Base.Six, ~CO2*Heat)

m_tukey=contrast(M1.Base.Six.post, method='pairwise', by=NULL)
summary(m_tukey)

####################--------------------------------------------------------------
#Shell Corrosion Data
####################--------------------------------------------------------------
setwd("/Users/racinerangel/Desktop/Climate Manipulations/Shell Analyses")
#----------------------------------------------------------------------------
#SHELL CORROSION ANALYSIS---------------------------------------------------
ShellInner<-read.csv("Shell_Corrosion_Inner.csv")
ShellOuter<-read.csv("Shell_Corrosion_Outer.csv")

#Transform outer shell to deal with negative skew
ShellOuter$SQPer_Corroded<-(sqrt(max(ShellOuter$Per_Corroded+1)-ShellOuter$Per_Corroded))

####################--------------------------------------------------------------
ShellInner$Treatment<-factor(ShellInner$Treatment, levels = c("Control", "Heat", "CO2", "Both"), labels=c("Unmanipulated", "Warming", "CO2 Addition", "Warming + CO2"))

ShellOuter$Treatment<-factor(ShellOuter$Treatment, levels = c("Control", "Heat", "CO2", "Both"), labels=c("Unmanipulated", "Warming", "CO2 Addition", "Warming + CO2"))


####################--------------------------------------------------------------
#Descriptive stats by Month and Pool Treatment
####################--------------------------------------------------------------
Corr.In %>% count(Treatment, Timepoint, Valve)
####################--------------------------------------------------------------
#Corrosion Averages
####################--------------------------------------------------------------
Corr.In <- ShellInner %>%
  group_by(Month,Timepoint, Treatment) %>%
  summarize(mean = mean(Per_Corroded),
            n = sum(!is.na(Treatment)),
            sd= sd(Per_Corroded),
            se=sd/sqrt(n))

Corr.Out <- ShellOuter %>%
  group_by(Month, Timepoint, Treatment) %>%
  summarize(mean = mean(Per_Corroded),
            n = sum(!is.na(Treatment)),
            sd= sd(Per_Corroded),
            se=sd/sqrt(n))

#Updated average files for figures
write.csv(Corr.In, "MS_AvgShellCorr_Inner_update.csv")
write.csv(Corr.Out, "MS_AvgShellCorr_Outer_update.csv")


####################--------------------------------------------------------------
#INNER SHELL CORROSION
####################--------------------------------------------------------------
Corr.Zero<-filter(ShellInner, Month == "March")
Corr.Three<-filter(ShellInner, Month =="July")
Corr.Six<-filter(ShellInner, Month =="September")

#Check distribution
hist(Corr.Zero$Gam_Corr) 
shapiro.test(Corr.Zero$Gam_Corr) #not normal p-value = 0.0005936 - Gamma

hist(Corr.Three$Gam_Corr)
shapiro.test(Corr.Three$Per_Corroded) # not normal  p-value =  0.0002151 - Gamma

hist(Corr.Six$Gam_Corr)
shapiro.test(Corr.Six$Gam_Corr) #not normal p-value = 0.004605 - Gamma

####################--------------------------------------------------------------
#GLM with Time 0
####################--------------------------------------------------------------
M.Corr.Zero<-glm(Gam_Corr~CO2*Heat+Volume, data=Corr.Zero, family=Gamma(link='log'))
summary(M.Corr.Zero)
anova(M.Corr.Zero, test="F")

#No difference
check_model(M.Corr.Zero) #decent

M.Corr.Zero.post<- emmeans(M.Corr.Zero, ~CO2*Heat)

m_tukey=contrast(M.Corr.Zero.post, method='pairwise', by=NULL)
summary(m_tukey) 


####################--------------------------------------------------------------
#GLM with Time 3
####################--------------------------------------------------------------
M.Corr.Three<-glm(Gam_Corr~CO2*Heat+Volume, data=Corr.Three, family=Gamma(link='log'))
summary(M.Corr.Three)
anova(M.Corr.Three, test="F")


check_model(M.Corr.Three) #decent

M.Corr.Three.post<- emmeans(M.Corr.Three, ~CO2*Heat)

m_tukey=contrast(M.Corr.Three.post, method='pairwise', by=NULL)
summary(m_tukey) #no SIG DIFFERENCES


####################--------------------------------------------------------------
#GLM with Time 6
####################--------------------------------------------------------------
M.Corr.Six<-glm(Gam_Corr~CO2*Heat+Volume, data=Corr.Six, family=Gamma(link='log'))
summary(M.Corr.Six)
anova(M.Corr.Six, test="F")

#CO2 Sig, CO2*Heat

check_model(M.Corr.Six) #decent

M.Corr.Six.post<- emmeans(M.Corr.Six, ~CO2*Heat)

m_tukey=contrast(M.Corr.Six.post, method='pairwise', by=NULL)
summary(m_tukey) #SIG DIFFERENCES


plot(Gam_Corr~Treatment, data=Corr.Six)

####################--------------------------------------------------------------
#OUTER SHELL CORROSION
####################--------------------------------------------------------------
Corr.Zero<-filter(ShellOuter, Month == "March")
Corr.Three<-filter(ShellOuter, Month =="July")
Corr.Six<-filter(ShellOuter, Month =="September")

#Check distribution
hist(Corr.Zero$SQPer_Corroded) 
shapiro.test(Corr.Zero$SQPer_Corroded) #normal  p-value = 0.4351

hist(Corr.Three$SQPer_Corroded)
shapiro.test(Corr.Three$Per_Corroded) #normal  p-value = 0.09117

hist(Corr.Six$SQPer_Corroded)
shapiro.test(Corr.Six$SQPer_Corroded) #normal p-value = 0.8795

hist(Corr.36$SQPer_Corroded)
shapiro.test(Corr.36$SQPer_Corroded) #Normal  p-value = 0.5203

hist(OuterCorr.avg$Mean.Corr)
shapiro.test(OuterCorr.avg$Mean.Corr) #p-value = 0.2796 normal

####################--------------------------------------------------------------
#GLM with Time 0 
####################--------------------------------------------------------------
M.Corr.Zero<-glm(SQPer_Corroded~CO2*Heat+Volume, data=Corr.Zero)
summary(M.Corr.Zero)
anova(M.Corr.Zero, test="F")
#No difference

check_model(M.Corr.Zero) #decent

M.Corr.Zero.post<- emmeans(M.Corr.Zero, ~CO2*Heat)

m_tukey=contrast(M.Corr.Zero.post, method='pairwise', by=NULL)
summary(m_tukey)

####################--------------------------------------------------------------
#GLM with Time 3 
####################--------------------------------------------------------------
M.Corr.Three<-glm(SQPer_Corroded~CO2*Heat+Volume, data=Corr.Three)
summary(M.Corr.Three)
anova(M.Corr.Three, test="F")
#No difference
Anova(M.Corr.Three, test="F", type=3)


check_model(M.Corr.Three) #decent

M.Corr.Three.post<- emmeans(M.Corr.Three, ~CO2*Heat)

m_tukey=contrast(M.Corr.Three.post, method='pairwise', by=NULL)
summary(m_tukey)

####################--------------------------------------------------------------
#GLM with Time 6 
####################--------------------------------------------------------------
M.Corr.Six <- glm(SQPer_Corroded~CO2*Heat+Volume, data=Corr.Six)
summary(M.Corr.Six)
anova(M.Corr.Six, test="F")
#CO2 Sig

check_model(M.Corr.Six) #decent

M.Corr.Six.post<- emmeans(M.Corr.Six, ~CO2*Heat)

m_tukey=contrast(M.Corr.Six.post, method='pairwise', by=NULL)
summary(m_tukey) #SIG DIFFERENCES


####################--------------------------------------------------------------
#Calculate biomass of pools
####################--------------------------------------------------------------
setwd("/Users/racinerangel/Desktop/Climate Manipulations")


Biomass<-read.csv("Mussel_Biomass_Jeb.csv")

Biomass$Treatment<-factor(Biomass$Treatment, levels = c("Control", "Heat", "CO2", "Both"), labels=c("Unmanipulated", "Warming", "CO2 Addition", "Warming + CO2"))

Biomass$Month<-factor(Biomass$Month, levels = c("March", "July", "September"))

Muss<- Biomass %>%
  group_by(Month, Treatment) %>%
  summarize(mean = mean(Norm_Biomass),
            sd = sd(Norm_Biomass))



