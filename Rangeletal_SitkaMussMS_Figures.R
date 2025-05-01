################################################################
#Title: Rangel et al. Sitka Mussel Climate Manipulation Figures
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
#############################################
#Set working Directory
#############################################
setwd("/Users/racinerangel/Desktop/Climate Manipulations")
#######################################################

####################--------------------------------------------------------------
#Environmental Data
####################--------------------------------------------------------------

####################--------------------------------------------------------------
#pH Data All Pools - Hanna Data
####################--------------------------------------------------------------
pH_Meta<- read.csv("pH_Meta.csv")
str(pH_Meta) #Check structure of data frame

Day<-filter(pH_Meta, DayNight == "Day" )
Night<-filter(pH_Meta, DayNight == "Night")     


####################--------------------------------------------------------------
#Day pH by MONTH - Figure 1
####################--------------------------------------------------------------
Day$Month<-factor(Day$Month, levels = c("March", "July", "September"), labels = c("0", "3", "6"))
Day$Treatment<-factor(Day$Treatment, levels = c("Unmanipulated", "CO2", "Warming", "Warming + CO2"), labels=c("Unmanipulated", "CO2 Added", "Warmed", "CO2 Added + Warmed"))


pH.day<-Day%>%
  ggplot(aes(x=Month, y=HannapH_Tot, fill=Treatment, col=Treatment)) + labs(y="Hanna pH") +
  geom_boxplot(position=position_dodge(0.85), pch=21, width = 0.75) + ylim(6,10) +  geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.85), pch=21, show.legend = F) +
  scale_fill_manual(name="Treatment", values=c("white","gray47", "white", "firebrick"), 
                    labels=c("Unmanipulated",
                             "CO2 Added" = bquote(CO[2]~"Added"),
                             "Warmed",
                             "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed"))) +
  scale_color_manual(name="Treatment", values=c("gray47", "black", "firebrick", "black"),
                     labels=c("Unmanipulated",
                              "CO2 Added" = bquote(CO[2]~"Added"),
                              "Warmed",
                              "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed")))+
  theme_bw() +
  theme(plot.margin =  margin(.5, .5, 0, .5, "cm"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.text.x  = element_text(size=16), axis.text.y  = element_text(size=16),
        axis.title  = element_text(size=16), legend.title  = element_text(size=16), legend.text=element_text(size=16)) +
  labs(x ="", y = "pH (Total Units)") + geom_vline(xintercept=1.5, linetype="dashed", color = "gray70")

pH.day

pH.day<- pH.day + ggtitle("Day") + theme(plot.title=element_text(hjust=0.5, vjust=0.5, size=22))


####################--------------------------------------------------------------
##Night pH by MONTH
####################--------------------------------------------------------------

Night$Month<-factor(Night$Month, levels = c("March", "July", "September"), labels = c("0", "3", "6"))
Night$Treatment<-factor(Night$Treatment, levels = c("Unmanipulated", "CO2", "Warming", "Warming + CO2"), labels=c("Unmanipulated", "CO2 Added", "Warmed", "CO2 Added + Warmed"))


pH.night<-Night%>%
  ggplot(aes(x=Month, y=HannapH_Tot, fill=Treatment, col=Treatment)) + labs(y="Hanna pH") +
  geom_boxplot(position=position_dodge(0.85), pch=21, show.legend = T, width = 0.75) + ylim(6,10) +  geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.85), pch=21, show.legend = F) +
  scale_fill_manual(name="Treatment", values=c("white","gray47", "white", "firebrick"), 
                    labels=c("Unmanipulated",
                             "CO2 Added" = bquote(CO[2]~"Added"),
                             "Warmed",
                             "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed"))) +
  scale_color_manual(name="Treatment", values=c("gray47", "black", "firebrick", "black"),
                     labels=c("Unmanipulated",
                              "CO2 Added" = bquote(CO[2]~"Added"),
                              "Warmed",
                              "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed")))+
  theme_bw() +
  theme(plot.margin =  margin(.5, .5, 0, .5, "cm"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.text.x  = element_text(size=16), axis.text.y  = element_text(size=16),
        axis.title  = element_text(size=16), legend.title  = element_text(size=16), legend.text=element_text(size=16)) +
  labs(x ="", y = "") + geom_vline(xintercept=1.5, linetype="dashed", color = "gray70")

pH.night 

pH.night <- pH.night + ggtitle("Night") + theme(plot.title=element_text(hjust=0.5, vjust=0.5, size=22))


####################--------------------------------------------------------------
#TEMPERATURE DATA
####################--------------------------------------------------------------
Temp<-read.csv("DetidedTemp.csv")

Temp$Month<-factor(Temp$Month, levels = c("0", "3", "6"), labels = c("0", "3", "6"))
Temp$Treatment<-factor(Temp$Treatment, levels = c("Unmanipulated", "CO2 Addition", "Warming", "Warming + CO2"), labels=c("Unmanipulated", "CO2 Added", "Warmed", "CO2 Added + Warmed"))

DayT<-filter(Temp, DayNight == "Day" )
NightT<-filter(Temp, DayNight == "Night")


####################--------------------------------------------------------------
#Temperature 99%
####################--------------------------------------------------------------

#99%
Temp.99<-DayT%>%
  ggplot(aes(x=Month, y=Daily99MaxPercentile, fill=Treatment, col=Treatment)) + labs(y="Temp") +
  geom_boxplot(position=position_dodge(0.85), pch=21, show.legend = T, width = 0.75) + ylim(5,25) +  geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.85), pch=21, show.legend = F) +
  scale_fill_manual(name="Treatment", values=c("white","gray47", "white", "firebrick"), 
                    labels=c("Unmanipulated",
                             "CO2 Added" = bquote(CO[2]~"Added"),
                             "Warmed",
                             "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed"))) +
  scale_color_manual(name="Treatment", values=c("gray47", "black", "firebrick", "black"),
                     labels=c("Unmanipulated",
                              "CO2 Added" = bquote(CO[2]~"Added"),
                              "Warmed",
                              "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed")))+
  theme_bw() +
  theme(plot.margin =  margin(0.3, .5, 0, .5, "cm"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.text.x  = element_text(size=16), axis.text.y  = element_text(size=16),
        axis.title.y  = element_text(size=14), legend.title  = element_text(size=16), legend.text=element_text(size=16)) +
  labs(x ="", y = "Daily 99th Percentile Pool Temperatures (°C)",) + geom_vline(xintercept=1.5, linetype="dashed", color = "gray70")

Temp.99

Temp.99 <- Temp.99 + ggtitle("Day") + theme(plot.title=element_text(hjust=0.5, vjust=0.5, size=22))



Temp.99Night<-NightT%>%
  ggplot(aes(x=Month, y=Daily99MaxPercentile, fill=Treatment, col=Treatment)) + labs(y="Temp") +
  geom_boxplot(position=position_dodge(0.85), pch=21, show.legend = T, width = 0.75) + ylim(5,25) +  geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.85), pch=21, show.legend = F) +
  scale_fill_manual(name="Treatment", values=c("white","gray47", "white", "firebrick"), 
                    labels=c("Unmanipulated",
                             "CO2 Added" = bquote(CO[2]~"Added"),
                             "Warmed",
                             "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed"))) +
  scale_color_manual(name="Treatment", values=c("gray47", "black", "firebrick", "black"),
                     labels=c("Unmanipulated",
                              "CO2 Added" = bquote(CO[2]~"Added"),
                              "Warmed",
                              "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed")))+
  theme_bw() +
  theme(plot.margin =  margin(0.3, .5, 0, .5, "cm"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.text.x  = element_text(size=16), axis.text.y  = element_text(size=16),
        axis.title.y  = element_text(size=14), legend.title  = element_text(size=16), legend.text=element_text(size=16)) +
  labs(x ="", y = "") + geom_vline(xintercept=1.5, linetype="dashed", color = "gray70")

Temp.99Night

Temp.99Night <- Temp.99Night + ggtitle("Night") + theme(plot.title=element_text(hjust=0.5, vjust=0.5, size=22))


#figure.DayNightTemp<-ggarrange(Temp.99, Temp.99Night,labels=("B"),
                               #common.legend = FALSE, legend = "none", font.label = list(size=20), ncol=2, nrow=1)
#figure.DayNightTemp<-annotate_figure(figure.DayNightTemp, bottom = text_grob("Time Of Exposure (Months)", size=16, hjust=0.8, vjust = -1.3))
#pdf("DayNightTemp.pdf",height=6.0,width=10)
#figure.DayNightTemp
#dev.off()
#Temp.99
#pdf("TempDaily99.pdf",height=5,width=7)
#Temp.99
#dev.off()


#Both pH and Temp together - Make sure pH Figures are made - Figure 1

figure.DayNight<-ggarrange(pH.day, pH.night,Temp.99, Temp.99Night, labels=c("A", "B", "C", "D"), vjust=1.3,hjust=-2.5,
                           common.legend = TRUE, legend = "right", font.label = list(size=20), ncol=2, nrow=2)


figure.DayNight<-annotate_figure(figure.DayNight, bottom = text_grob("Time Of Exposure (Months)", size=16, hjust=0.8, vjust = -1.0))


ggsave("TempandpHDayNightV2.svg", figure.DayNight, device="svg")
ggsave("TempandpHDayNightV2.pdf",height=18, width=22, units = c("cm"), dpi=1200, useDingbats=FALSE)

figure.DayNight 

dev.off()



####################--------------------------------------------------------------
#SHELL THICKNESS - Figure 2
####################--------------------------------------------------------------
setwd("/Users/racinerangel/Desktop/Climate Manipulations/Shell Analyses/Shell Thickness")
####################--------------------------------------------------------------
ShellTreatA<-read.csv("avgATreatment.csv") #Only used for plotting
ShellTreatB<-read.csv("avgBTreatment.csv")
ShellTreatC<-read.csv("avgCtreatment.csv")

ShellTreatA$Timepoint<-factor(ShellTreatA$Timepoint, levels = c("0", "3", "6"), labels = c("0", "3", "6"))
ShellTreatB$Timepoint<-factor(ShellTreatB$Timepoint, levels = c("0", "3", "6"), labels = c("0", "3", "6"))
ShellTreatC$Timepoint<-factor(ShellTreatC$Timepoint, levels = c("0", "3", "6"), labels = c("0", "3", "6"))

ShellTreatA$Treatment<-factor(ShellTreatA$Treatment, levels = c("Control", "CO2", "Heat", "Both"), labels=c("Unmanipulated", "CO2 Added","Warmed", "CO2 Added + Warmed")) 

ShellTreatB$Treatment<-factor(ShellTreatB$Treatment, levels = c("Control", "CO2", "Heat", "Both"), labels=c("Unmanipulated", "CO2 Added", "Warmed", "CO2 Added + Warmed"))

ShellTreatC$Treatment<-factor(ShellTreatC$Treatment, levels = c("Control", "CO2", "Heat", "Both"), labels=c("Unmanipulated", "CO2 Added","Warmed", "CO2 Added + Warmed")) 

####################--------------------------------------------------------------
#LIP SHELL THICKNESS
####################--------------------------------------------------------------
#jpeg("lineslipSE.jpg",height=550*4,width=650*4,res=100*4)

Lip<-ggplot(ShellTreatA, aes(x=Timepoint, y=mean, group=Treatment)) +
  geom_point(aes(shape=Treatment, color=Treatment), size=4.5, stroke=1.8, position=position_dodge(width=0.2))+
  scale_shape_manual(values=c(1, 16, 1, 16), labels=c("Unmanipulated",
                                                      "CO2 Added" = bquote(CO[2]~"Added"),
                                                      "Warmed",
                                                      "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed")))+  
  scale_color_manual(values=c("gray47", "gray47","firebrick", "firebrick"), labels=c("Unmanipulated",
                                                                                     "CO2 Added" = bquote(CO[2]~"Added"),
                                                                                     "Warmed",
                                                                                     "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed"))) + geom_line(linetype="solid", position=position_dodge(width=0.2))+theme_bw()+
  theme(plot.margin =  margin(.5, 0, .5, .5, "cm"), legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.text.x  = element_text(size=16), axis.text.y  = element_text(size=16),
        axis.title  = element_text(size=16), legend.title  = element_text(size=16), legend.text=element_text(size=16), plot.title=element_text(size=16, hjust=0.5)) +
  labs(x ="", y = "Average Shell Thickness (mm)", title="Lip")+ geom_errorbar(aes(ymin=mean-sd, 
                                                              ymax=mean+sd), alpha=0.6, width=.5, position=position_dodge(0.2)) + scale_y_continuous(name="Standardized Shell Thickness (mm)", limits=c(0.000, 0.050))+
  annotate("text", x=c(1,2,3), y=c(0.041, 0.041, 0.041), label=c("", "",""), size=12) + geom_vline(xintercept=1.5, linetype="dashed", color = "gray70")

Lip

dev.off()  


####################--------------------------------------------------------------
#MIDDLE SHELL THICKNESS
####################--------------------------------------------------------------
#jpeg("linesmidSE.jpg",height=550*4,width=650*4,res=100*4)

Mid<-ggplot(ShellTreatB, aes(x=Timepoint, y=mean, group=Treatment)) +
  geom_point(aes(shape=Treatment, color=Treatment), size=4.5, stroke=1.8, position=position_dodge(0.2))+
  scale_shape_manual(values=c(1, 16, 1, 16), labels=c("Unmanipulated",
                                                      "CO2 Added" = bquote(CO[2]~"Added"),
                                                      "Warmed",
                                                      "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed")))+  
  scale_color_manual(values=c("gray47", "gray47","firebrick", "firebrick"), labels=c("Unmanipulated",
                                                                                     "CO2 Added" = bquote(CO[2]~"Added"),
                                                                                     "Warmed",
                                                                                     "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed")))+
  geom_line(linetype="solid", position=position_dodge(0.2)) + theme_bw()+
  theme(plot.margin =  margin(.5, .5, .5, 0, "cm"),axis.line = element_line(colour = "black"), legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.text.x  = element_text(size=16), axis.text.y  = element_text(size=16),
        axis.title  = element_text(size=16), legend.title  = element_text(size=16), legend.text=element_text(size=16), plot.title=element_text(size=16, hjust=0.5))+
  labs(x ="", y = "", title="Middle") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), alpha=0.6, width=.5, position=position_dodge(0.2))+ 
  scale_y_continuous(name="", limits=c(0.000, 0.043))+annotate("text", x=c(1,2,3), y=c(0.041, 0.041, 0.041), label=c("", ""," "), size=12) + geom_vline(xintercept=1.5, linetype="dashed", color = "gray70")

Mid

dev.off()

####################--------------------------------------------------------------
#BASE SHELL THICKNESS
####################--------------------------------------------------------------
#jpeg("linesbaseSE.jpg",height=550*4,width=650*4,res=100*4)

Base<-ggplot(ShellTreatC, aes(x=Timepoint, y=mean, group=Treatment)) +
  geom_point(aes(shape=Treatment, color=Treatment), size=4.5, stroke=1.8, position=position_dodge(0.2))+
  scale_shape_manual(values=c(1, 16, 1, 16), labels=c("Unmanipulated",
                                                      "CO2 Added" = bquote(CO[2]~"Added"),
                                                      "Warmed",
                                                      "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed")))+  
  scale_color_manual(values=c("gray47", "gray47","firebrick", "firebrick"), labels=c("Unmanipulated",
                                                                                     "CO2 Added" = bquote(CO[2]~"Added"),
                                                                                     "Warmed",
                                                                                      "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed")))+
  geom_line(linetype="solid", position=position_dodge(0.2)) + theme_bw()+
  theme(plot.margin =  margin(.5, .5, .5, 0, "cm"),axis.line = element_line(colour = "black"), legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.text.x  = element_text(size=16), axis.text.y  = element_text(size=16),
        axis.title  = element_text(size=16), legend.title  = element_text(size=16), legend.text=element_text(size=16), plot.title=element_text(size=16, hjust=0.5))+
  labs(x ="", y = "", title="Base") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), alpha=0.6, width=.5, position=position_dodge(0.2))+ 
  scale_y_continuous(name="", limits=c(0.000, 0.043))+annotate("text", x=c(1,2,3), y=c(0.041, 0.041, 0.041), label=c("", ""," "), size=12) + geom_vline(xintercept=1.5, linetype="dashed", color = "gray70")

Base

dev.off()

####################--------------------------------------------------------------
#Full Shell Thickness Plot - Figure 2
####################--------------------------------------------------------------
figure.Thickness<-ggarrange(Lip, Mid, Base, labels= c("A", "B", "C"),
                            common.legend = TRUE, legend = "right", font.label = list(size=20), ncol=3, nrow=1)

figure.Thickness<-annotate_figure(figure.Thickness, bottom = text_grob("Time Of Exposure (Months)", size=16, hjust=0.8, vjust = -1.3))

quartz()
#pdf("Thickness.pdf",height=5,width=14)
#jpeg("Thickness.jpg",height=550*4,width=1250*4,res=100*5)


figure.Thickness

dev.off()

#Have to use Cairo for PDF to show graphics
ggsave("thickness.pdf", figure.Thickness, width=14, height=5, dpi=900, device = cairo_pdf)
ggsave("thickness.svg", figure.Thickness, device="svg")



####################--------------------------------------------------------------
#Shell Strength Data - Figure 3
####################--------------------------------------------------------------
setwd("/Users/racinerangel/Desktop/Climate Manipulations/Shell analyses")
####################--------------------------------------------------------------
Muss_Str <- read.csv("ShellStrength_CL.csv")

Muss_Str<- Muss_Str %>% filter(Pool!="21") #Removing Pool 21 from analyses

#DO NOT FORGET***************Reorder treatments
Muss_Str$Treatment<-factor(Muss_Str$Treatment, levels = c("Control","CO2","Heat", "Both"), labels = c("Unmanipulated", "CO2 Added", "Warmed", "CO2 Added + Warmed"))

Muss_Str$Month<-factor(Muss_Str$Month, levels = c("March", "July", "September"), labels = c("March", "July", "September"))


####################--------------------------------------------------------------
#SHELL STRENGTH - FIGURE 
####################--------------------------------------------------------------
Str<-ggplot(Strength.avg, aes(x=Month, y=mean, group=Treatment)) +
  geom_point(aes(shape=Treatment, color=Treatment), size=4.5, stroke=1.8, position=position_dodge(width=0.2))+
  scale_shape_manual(values=c(1, 16, 1, 16), labels=c("Unmanipulated",
                                                      "CO2 Added" = bquote(CO[2]~"Added"),
                                                      "Warmed",
                                                      "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed")))+  
  scale_color_manual(values=c("gray47", "gray47","firebrick", "firebrick"), labels=c("Unmanipulated",
                                                                                     "CO2 Added" = bquote(CO[2]~"Added"),
                                                                                     "Warmed",
                                                                                     "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed"))) + geom_line(linetype="solid", position=position_dodge(width=0.2))+theme_bw()+
  theme(plot.margin =  margin(0.5, 0.5, 0.5, 0.5, "cm"), legend.position = "right",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.text.x  = element_text(size=16), axis.text.y  = element_text(size=16),
        axis.title  = element_text(size=16), legend.title  = element_text(size=16), legend.text=element_text(size=16), plot.title=element_text(size=16, hjust=0.5)) +
  labs(x ="Time Of Exposure (Months)", y = "Average Standardized Shell Strength (N)")+ geom_errorbar(aes(ymin=mean-sd, 
                                                                                  ymax=mean+sd), alpha=0.6, width=.2, position=position_dodge(0.2)) + scale_y_continuous(name="Average Standardized Shell Strength (N)", limits=c(0.000, 3.50)) + geom_vline(xintercept=1.5, linetype="dashed", color = "gray70") + scale_x_discrete(breaks=c("March","July","September"), labels=c("0","3", "6"))

Str

ggsave("Strength_Line.svg", Str, device="svg")
ggsave("Strength.Line.pdf",height=16, width=26, units = c("cm"), dpi=1200, useDingbats=FALSE)



####################--------------------------------------------------------------
#Shell Corrosion Data
####################--------------------------------------------------------------
setwd("/Users/racinerangel/Desktop/Climate Manipulations/Shell Analyses")
#----------------------------------------------------------------------------
#SHELL CORROSION ANALYSIS---------------------------------------------------
ShellInner<-read.csv("Shell_Corrosion_Inner.csv")
ShellOuter<-read.csv("Shell_Corrosion_Outer.csv")
InnerCorr.avg<- read.csv("MS_AvgShellCorr_Inner_update.csv")
OuterCorr.avg<- read.csv("MS_AvgShellCorr_Outer_update.csv")
####################--------------------------------------------------------------
InnerCorr.avg$Treatment<-factor(InnerCorr.avg$Treatment, levels = c("Unmanipulated", "CO2 Added","Warmed", "CO2 Added + Warmed"), labels=c("Unmanipulated", "CO2 Added","Warmed", "CO2 Added + Warmed"))

OuterCorr.avg$Treatment<-factor(OuterCorr.avg$Treatment, levels = c("Control", "CO2","Heat", "Both"), labels=c("Unmanipulated", "CO2 Added", "Warmed", "CO2 Added + Warmed"))


####################--------------------------------------------------------------
#INNER SHELL CORROSION GEOMPOINT
####################--------------------------------------------------------------
Inner<-InnerCorr.avg%>%
  mutate(Month= fct_relevel(Month, "March", "July","September"))%>%
  ggplot(aes(x=Month, y=mean, group=Treatment)) +
  geom_point(aes(shape=Treatment, color=Treatment), size=4.5, stroke=1.8, position=position_dodge(0.2))+
  scale_shape_manual(values=c(1, 16, 1, 16), labels=c("Unmanipulated",
                                                      "CO2 Added" = bquote(CO[2]~"Added"),
                                                      "Warmed",
                                                      "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed")))+  
  scale_color_manual(values=c("gray47", "gray47","firebrick", "firebrick"), labels=c("Unmanipulated",
                                                                                     "CO2 Added" = bquote(CO[2]~"Added"),
                                                                                     "Warmed",
                                                                                     "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed")))+
  geom_line(linetype="solid", position=position_dodge(0.2)) + theme_bw()+
  theme(plot.margin =  margin(.5, .5, .5, 0.5, "cm"),axis.line = element_line(colour = "black"), legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.text.x  = element_text(size=16), axis.text.y  = element_text(size=16),
        axis.title  = element_text(size=16), legend.title  = element_text(size=16), legend.text=element_text(size=16), plot.title=element_text(size=16, hjust=0.5))+
  labs(x ="", y = "Shell Corrosion (%)", title="") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), alpha=0.6, width=.5, position=position_dodge(0.2))+ geom_vline(xintercept=1.5, linetype="dashed", color = "gray70")+ scale_x_discrete(breaks=c("March","July","September"), labels=c("0","3", "6")) +scale_y_continuous(name="Shell Corrosion (%)", limits = c(-15, 120), breaks = seq(0, 120, by = 20))


Inner


####################--------------------------------------------------------------
#OUTERSHELL CORROSION GEOMPOINT
####################--------------------------------------------------------------


Outer<-OuterCorr.avg%>%
  mutate(Month= fct_relevel(Month, "March", "July","September"))%>%
  ggplot(aes(x=Month, y=mean, group=Treatment)) +
  geom_point(aes(shape=Treatment, color=Treatment), size=4.5, stroke=1.8, position=position_dodge(0.2))+
  scale_shape_manual(values=c(1, 16, 1, 16), labels=c("Unmanipulated",
                                                      "CO2 Added" = bquote(CO[2]~"Added"),
                                                      "Warmed",
                                                      "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed")))+  
  scale_color_manual(values=c("gray47", "gray47","firebrick", "firebrick"), labels=c("Unmanipulated",
                                                                                     "CO2 Added" = bquote(CO[2]~"Added"),
                                                                                     "Warmed",
                                                                                     "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed")))+
  geom_line(linetype="solid", position=position_dodge(0.2)) + theme_bw()+
  theme(plot.margin =  margin(.5, .5, 0.5, 0.5, "cm"),axis.line = element_line(colour = "black"), legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.text.x  = element_text(size=16), axis.text.y  = element_text(size=16),
        axis.title  = element_text(size=16), legend.title  = element_text(size=16), legend.text=element_text(size=16), plot.title=element_text(size=16, hjust=0.5))+
  labs(x ="", y = "", title="") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), alpha=0.6, width=.5, position=position_dodge(0.2))+ geom_vline(xintercept=1.5, linetype="dashed", color = "gray70")+ scale_x_discrete(breaks=c("March","July","September"), labels=c("0","3", "6")) +scale_y_continuous(name="Shell Corrosion (%)", limits = c(-15, 120), breaks = seq(0, 120, by = 20))

Outer

figure.InOut<-ggarrange(Outer, Inner, labels = c("A", "B"),
                        common.legend = TRUE, legend = "right", font.label = list(size=20), ncol=1, nrow=2)


figure.InOut<-annotate_figure(figure.InOut, bottom = text_grob("Time Of Exposure (Months)", size=16, hjust=0.8, vjust = -1.3))


pdf("ShellCorrosion.pdf",height=14,width=10,)
ggsave("ShellCorrosionline.pdf",height=14, width=10, dpi=900, device = cairo_pdf)
ggsave("ShellCorrosionline.svg", figure.InOut, device="svg")

figure.InOut

dev.off()


####################--------------------------------------------------------------
#SUPPLEMENTARY FIGURES
####################--------------------------------------------------------------

####################--------------------------------------------------------------
#Metabolic Rate Anaylsis
####################--------------------------------------------------------------
setwd("/Users/racinerangel/Desktop/Climate Manipulations/Climate Manipulations MR")
####################--------------------------------------------------------------
Q10_Meta<-read.csv("Q10_M_MJS_CL.csv")
Q10.avg<- read.csv("Q10.Avg.csv") #Averaged values

Q10_Meta<- Q10_Meta %>% filter(PoolID!="Pool21")

####################--------------------------------------------------------------
####################--------------------------------------------------------------
Q10_Meta_Three<-filter(Q10_Meta, ThreeTrials == "Y")
count(Q10_Meta$Treatment)

#Q10_Meta$Month<-factor(Q10_Meta$Month, levels = c("March", "July", "Sept"), labels = c("0", "3", "6")) #For plotting


#DO NOT FORGET***************Reorder treatments
Q10_Meta$Treatment<-factor(Q10_Meta$Treatment, levels = c("Control", "Heat", "CO2", "Both"), labels = c("Unmanipulated", "Warming", "CO2 Addition", "Warming + CO2"))

####################--------------------------------------------------------------
#Q10 Plot
####################--------------------------------------------------------------

QTime<-Q10_Meta%>%
  mutate(Month = fct_relevel(Month, "March", "July","Sept"))%>%
  ggplot(aes(x=Month, y=Q10, fill=Treatment, col=Treatment)) + labs(y=expression(Q[10])) +
  geom_boxplot(position=position_dodge(1)) + geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 1), pch=21, show.legend = F) +
  scale_fill_manual(name="Treatment", values=c("white", "white","gray47","firebrick"), 
                    labels=c("Unmanipulated",
                             "Warming",
                             "CO2 Addition" = bquote(CO[2]~"Addition"),
                             "Warming + CO2" = bquote("Warming + CO" [2]))) +
  scale_color_manual(name="Treatment", values=c("gray47", "firebrick","black","black"),
                     labels=c("Unmanipulated",
                              "Warming",
                              "CO2 Addition" = bquote(CO[2]~"Addition"),
                              "Warming + CO2" = bquote("Warming + CO" [2])))+
  theme_bw() +
  theme(plot.margin =  margin(.5, .5, .5, .5, "cm"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.text.x  = element_text(size=16), axis.text.y  = element_text(size=16),
        axis.title  = element_text(size=16), legend.title  = element_text(size=16), legend.text=element_text(size=16)) + scale_x_discrete(breaks=c("March","July","Sept"), labels=c("0","3","6")) +
  labs(x ="Time Of Exposure (Months)", y = (expression(Thermal~Sensitivity~italic(Q[10])))) + ylim(0,2.0) + geom_vline(xintercept=1.51, linetype="dashed", color = "gray42")

QTime

pdf("Q10.pdf",height=5,width=7)
QTime

dev.off()

####################--------------------------------------------------------------
#Day Calcite by MONTH
####################--------------------------------------------------------------

Cal.day<-Day%>%
  ggplot(aes(x=Month, y=Calcite, fill=Treatment, col=Treatment)) + labs(y="Calcite") +
  geom_boxplot(position=position_dodge(0.85), pch=21, width = 0.75) + ylim(-1,12) +  geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.85), pch=21, show.legend = F) +
  scale_fill_manual(name="Treatment", values=c("white","gray47", "white", "firebrick"), 
                    labels=c("Unmanipulated",
                             "CO2 Added" = bquote(CO[2]~"Added"),
                             "Warmed",
                             "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed"))) +
  scale_color_manual(name="Treatment", values=c("gray47", "black", "firebrick", "black"),
                     labels=c("Unmanipulated",
                              "CO2 Added" = bquote(CO[2]~"Added"),
                              "Warmed",
                              "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed")))+
  theme_bw() +
  theme(plot.margin =  margin(.5, .5, 0, .5, "cm"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.text.x  = element_text(size=16), axis.text.y  = element_text(size=16),
        axis.title  = element_text(size=16), legend.title  = element_text(size=16), legend.text=element_text(size=16)) +
  labs(x ="", y = "Calicite") + geom_vline(xintercept=1.5, linetype="dashed", color = "gray70")

Cal.day

Cal.day<- Cal.day + ggtitle("Day") + theme(plot.title=element_text(hjust=0.5, vjust=0.5, size=24))

pdf("CalDay.pdf",height=5,width=7)

dev.off()

####################--------------------------------------------------------------
##Night Calcite by MONTH
####################--------------------------------------------------------------

Night$Month<-factor(Night$Month, levels = c("March", "July", "September"), labels = c("0", "3", "6"))
Night$Treatment<-factor(Night$Treatment, levels = c("Unmanipulated", "Warming", "CO2", "Warming + CO2"), labels=c("Unmanipulated", "Warming", "CO2 Addition", "Warming + CO2"))


Cal.night<-Night%>%
  ggplot(aes(x=Month, y=Calcite, fill=Treatment, col=Treatment)) + labs(y="Calcite") +
  geom_boxplot(position=position_dodge(0.85), pch=21, show.legend = T, width = 0.75) + ylim(-0.7,11) +  geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.85), pch=21, show.legend = F) +
  scale_fill_manual(name="Treatment", values=c("white","gray47", "white", "firebrick"), 
                    labels=c("Unmanipulated",
                             "CO2 Added" = bquote(CO[2]~"Added"),
                             "Warmed",
                             "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed"))) +
  scale_color_manual(name="Treatment", values=c("gray47", "black", "firebrick", "black"),
                     labels=c("Unmanipulated",
                              "CO2 Added" = bquote(CO[2]~"Added"),
                              "Warmed",
                              "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed")))+
  theme_bw() +
  theme(plot.margin =  margin(.5, .5, 0, .5, "cm"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.text.x  = element_text(size=16), axis.text.y  = element_text(size=16),
        axis.title  = element_text(size=16), legend.title  = element_text(size=16), legend.text=element_text(size=16)) +
  labs(x ="", y = "") + geom_vline(xintercept=1.5, linetype="dashed", color = "gray70") +geom_hline(yintercept=0, linetype="solid", color = "gray70")

Cal.night 

Cal.night <- Cal.night + ggtitle("Night") + theme(plot.title=element_text(hjust=0.5, vjust=0.5, size=24))



#Both pH and Temp together
figure.DayNight.Cal<-ggarrange(Cal.day, Cal.night,
                               common.legend = TRUE, legend = "right", font.label = list(size=20), ncol=2, nrow=1)


figure.DayNight.Cal<-annotate_figure(figure.DayNight.Cal, bottom = text_grob("Time Of Exposure (Months)", size=16, hjust=0.8, vjust = -1.3))

pdf("TempandpHDayNight.pdf",height=10,width=11)

figure.DayNight.Cal

dev.off()


####################--------------------------------------------------------------
#Day Aragonite by MONTH
####################--------------------------------------------------------------

Arg.day<-Day%>%
  ggplot(aes(x=Month, y=Aragonite, fill=Treatment, col=Treatment)) + labs(y="Aragonite") +
  geom_boxplot(position=position_dodge(0.85), pch=21, width = 0.75) + ylim(-1,12) +  geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.85), pch=21, show.legend = F) +
  scale_fill_manual(name="Treatment", values=c("white","gray47", "white", "firebrick"), 
                    labels=c("Unmanipulated",
                             "CO2 Added" = bquote(CO[2]~"Added"),
                             "Warmed",
                             "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed"))) +
  scale_color_manual(name="Treatment", values=c("gray47", "black", "firebrick", "black"),
                     labels=c("Unmanipulated",
                              "CO2 Added" = bquote(CO[2]~"Added"),
                              "Warmed",
                              "CO2 Added + Warmed" = bquote(CO[2]~"Added + Warmed")))+
  theme_bw() +
  theme(plot.margin =  margin(.5, .5, 0, .5, "cm"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.text.x  = element_text(size=16), axis.text.y  = element_text(size=16),
        axis.title  = element_text(size=16), legend.title  = element_text(size=16), legend.text=element_text(size=16)) +
  labs(x ="Time Of Exposure (Months)", y = "Aragonite") + geom_vline(xintercept=1.5, linetype="dashed", color = "gray70")

Arg.day

Arg.day<- Arg.day + ggtitle("Day") + theme(plot.title=element_text(hjust=0.5, vjust=0.5, size=24))

pdf("ArgDay.pdf",height=5,width=7)

dev.off()

####################--------------------------------------------------------------
##Night Aragonite by MONTH
####################--------------------------------------------------------------

Night$Month<-factor(Night$Month, levels = c("March", "July", "September"), labels = c("0", "3", "6"))
Night$Treatment<-factor(Night$Treatment, levels = c("Unmanipulated", "Warming", "CO2", "Warming + CO2"), labels=c("Unmanipulated", "Warming", "CO2 Addition", "Warming + CO2"))


Arg.night<-Night%>%
  ggplot(aes(x=Month, y=Aragonite, fill=Treatment, col=Treatment)) + labs(y="Aragonite") +
  geom_boxplot(position=position_dodge(0.85), pch=21, show.legend = T, width = 0.75) + ylim(-1, 7.5) +  geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.85), pch=21, show.legend = F) +
  scale_fill_manual(name="Treatment", values=c("white", "white","gray47","firebrick"), 
                    labels=c("Unmanipulated",
                             "Warming",
                             "CO2 Addition" = bquote(CO[2]~"Addition"),
                             "Warming + CO2" = bquote("Warming + CO" [2]))) +
  scale_color_manual(name="Treatment", values=c("gray47", "firebrick","black","black"),
                     labels=c("Unmanipulated",
                              "Warming",
                              "CO2 Addition" = bquote(CO[2]~"Addition"),
                              "Warming + CO2" = bquote("Warming + CO" [2])))+
  theme_bw() +
  theme(plot.margin =  margin(.5, .5, 0, .5, "cm"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.text.x  = element_text(size=16), axis.text.y  = element_text(size=16),
        axis.title  = element_text(size=16), legend.title  = element_text(size=16), legend.text=element_text(size=16)) +
  labs(x ="", y = "") + geom_vline(xintercept=1.5, linetype="dashed", color = "gray70") +geom_hline(yintercept=0, linetype="solid", color = "gray70")

Arg.night 

Arg.night <- Arg.night + ggtitle("Night") + theme(plot.title=element_text(hjust=0.5, vjust=0.5, size=24))



#Both pH and Temp together
figure.DayNight.Calc<-ggarrange(Cal.day, Cal.night,Arg.day, Arg.night, labels=c("A", "", "B", ""),
                                common.legend = TRUE, legend = "right", font.label = list(size=20), ncol=2, nrow=2)

figure.DayNight.Calc<-annotate_figure(figure.DayNight.Calc, bottom = text_grob("Time Of Exposure (Months)", size=16, hjust=0.8, vjust = -1.3))

pdf("ArgandCalDayNight.pdf",height=10,width=11)

figure.DayNight.Calc

dev.off()

####################--------------------------------------------------------------
#Temperature Range
####################--------------------------------------------------------------

Temp.Range<-Temp%>%
  ggplot(aes(x=Month, y=DailyrangeTemp, fill=Treatment, col=Treatment)) + labs(y="Daily Range of Pool Temperatures (°C)") +
  geom_boxplot(position=position_dodge(0.85), pch=21, show.legend = T, width = 0.75) + ylim(0,13) + geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.85), pch=21, show.legend = F)+
  scale_fill_manual(name="Treatment", values=c("white", "white","gray47","firebrick"), 
                    labels=c("Unmanipulated",
                             "Warming",
                             "CO2 Addition" = bquote(CO[2]~"Addition"),
                             "Warming + CO2" = bquote("Warming + CO" [2]))) +
  scale_color_manual(name="Treatment", values=c("gray47", "firebrick","black","black"),
                     labels=c("Unmanipulated",
                              "Warming",
                              "CO2 Addition" = bquote(CO[2]~"Addition"),
                              "Warming + CO2" = bquote("Warming + CO" [2])))+
  theme_bw() +
  theme(plot.margin =  margin(.5, 0, .5, .5, "cm"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.text.x  = element_text(size=16), axis.text.y  = element_text(size=16),
        axis.title  = element_text(size=16), legend.title  = element_text(size=16), legend.text=element_text(size=16)) +
  labs(x ="Time Of Exposure (Months)", y = "Daily Range of Pool Temperatures (°C)") + geom_vline(xintercept=1.5, linetype="dashed", color = "gray42")

Temp.Range

pdf("TempDailyRange.pdf",height=5,width=7)

Temp.Range

dev.off()

####################--------------------------------------------------------------
#Temperature Maximum
####################--------------------------------------------------------------
Temp.Max<-Temp%>%
  ggplot(aes(x=Month, y=maxTemp, fill=Treatment, col=Treatment)) + labs(y="Maximum of Pool Temperatures (°C)") +
  geom_boxplot(position=position_dodge(0.85), pch=21, show.legend = T, width = 0.75) + ylim(5,40) + geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.85), pch=21, show.legend = F) +
  scale_fill_manual(name="Treatment", values=c("white", "white","gray47","firebrick"), 
                    labels=c("Unmanipulated",
                             "Warming",
                             "CO2 Addition" = bquote(CO[2]~"Addition"),
                             "Warming + CO2" = bquote("Warming + CO" [2]))) +
  scale_color_manual(name="Treatment", values=c("gray47", "firebrick","black","black"),
                     labels=c("Unmanipulated",
                              "Warming",
                              "CO2 Addition" = bquote(CO[2]~"Addition"),
                              "Warming + CO2" = bquote("Warming + CO" [2])))+
  theme_bw() +
  theme(plot.margin =  margin(.5, 0, .5, .5, "cm"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.text.x  = element_text(size=16), axis.text.y  = element_text(size=16),
        axis.title  = element_text(size=16), legend.title  = element_text(size=16), legend.text=element_text(size=16)) +
  labs(x ="Time Of Exposure (Months)", y = "Maximum of Pool Temperatures (°C)") + geom_vline(xintercept=1.5, linetype="dashed", color = "gray42")

Temp.Max

pdf("TempMaximum.pdf",height=5,width=7)

Temp.Range

dev.off()



####################--------------------------------------------------------------
#Community Data Bar Plot
####################--------------------------------------------------------------
setwd("/Users/racinerangel/Desktop/Climate Manipulations")
####################--------------------------------------------------------------

Community_Meta<-read.csv("Community_Meta_MonthMusselsTogether.csv")

#Full_Comm_Month<-Community_Meta %>% filter(Month!="April")
Full_Comm_Month<-Community_Meta %>% filter(Month!="May")
Full_Comm_Month<-Full_Comm_Month %>% filter(Month!="June")
Full_Comm_Month<-Full_Comm_Month %>% filter(Month!="August")

Full_Comm_Month$Month<-factor(Full_Comm_Month$Month, levels = c("March", "July", "September"), labels= c("0","3", "6"))

####################--------------------------------------------------------------
#Community w bar plot by month and treatment
####################--------------------------------------------------------------

#Full_Comm_Month$Month<-factor(Full_Comm_Month$Month, levels = c("March", "May", "June", "July", "August", "September"))

Full_Comm_Month$Treatment<-factor(Full_Comm_Month$Treatment, levels = c("Unmanipulated", "CO2 Added", "Warmed", "CO2 Added + Warmed"))

Full_Comm_Month$Species<-factor(Full_Comm_Month$Species, levels = c("Bare Rock","Barnacles","Dead Barnacles", "Mussel","Dead Mussels", "Coralline Algae","Crustose Coralline", "Other Sessile Calcifying","Other Sessile Non-Calcifying", "Brown Algae","Green Algae","Red Algae"),
                                labels =c("Bare Rock", "Barnacles", "Dead Barnacles","Mussels","Dead Mussels","Coralline Algae","Crustose Algae", "Other Sessile Calcifying","Other Sessile Non-Calcifying", 
                                          "Brown Algae","Green Algae","Red Algae")) 
Full_Comm_Month$Treatment = factor(Full_Comm_Month$Treatment, labels=c(
  "Unmanipulated",
  "CO[2]~Added",
  "Warmed",
  "CO[2]~Added + Warmed"))

#Figure by month----------------------------------------------
Comm.bar_Month<- ggplot(Full_Comm_Month, aes(fill = Species, x=Month, y = avg_percent_cover)) + ylim(0,135)+
  geom_col(position="stack", width=0.7) + 
  theme(axis.text.x = element_text(angle = 0, size = 14, color = "black", hjust=0.5, face= "bold"), axis.title.x = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 16, face = "bold"), legend.title = element_text(size = 16, face = "bold"), 
        legend.text = element_text(size = 12, face="bold", colour = "black"), legend.text.align = 0,
        axis.text.y = element_text(colour = "black", size = 12, face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(2, 0.5, 0.5, 2, "cm")) + labs(x="Month", y="Average Percent Cover") +
  scale_fill_manual(values=c("gray", "#0F8299", "#7ABECC", "#3E9FB3", "#B8DEE6", "#CC7A88", "#B33E52","#99600F", "#B3823E","#CCAA7A", "#A3CC7A", "#990F26")) + facet_grid(~ Treatment, labeller = label_parsed) + theme(strip.text = element_text(size=16))                                                                                                                      

Comm.bar_Month

quartz()

#png("CommunityBarPlotTreatmentChp3.png",height=550*4,width=650*4,res=100*4)
#pdf("CommunityBarPlot.pdf",height=6,width=7)
#Save as scg to easily edit in Power point
ggsave("Month_Community_MusselsTogether.svg", Comm.bar_Month, device="svg")

Comm.bar

dev.off()

####################--------------------------------------------------------------
#Temperature data across time
####################--------------------------------------------------------------
setwd("/Users/racinerangel/Desktop/Temperature")
####################--------------------------------------------------------------

#### TEMPERATURE FIGURE- Supplementary Figure 5 with only mussel pools -------------------------------------------------------------------------
setwd("/Users/racinerangel/Desktop/Temperature")
#DayTotalT<-read.csv("FullDaytimeLowTempCL.csv") #Use for detided Temp

FullT<-read.csv("FullTempCL_Updated.csv", stringsAsFactors = FALSE)

#Removing non mussel pools
FullT<-FullT[,!names(FullT) %in% c("Pool.20", "Pool.25", "Pool.26", "Pool.29", "Pool.30","Pool.21")]

#Check Lowest and highest value for figures
min(FullT$Pool.5, na.rm = TRUE)
min(FullT$Pool.7, na.rm = TRUE)
min(FullT$Pool.8, na.rm = TRUE)
min(FullT$Pool.9, na.rm = TRUE)
min(FullT$Pool.10, na.rm = TRUE)
min(FullT$Pool.11, na.rm = TRUE)
min(FullT$Pool.13, na.rm = TRUE)
min(FullT$Pool.15, na.rm = TRUE) #-1.015
min(FullT$Pool.16, na.rm = TRUE)
min(FullT$Pool.27, na.rm = TRUE)
min(FullT$Pool.31, na.rm = TRUE)
min(FullT$Pool.33, na.rm = TRUE)
min(FullT$Pool.35, na.rm = TRUE) #-1.356
min(FullT$Pool.36, na.rm = TRUE)

#Check Lowest and highest value for figures
max(FullT$Pool.5, na.rm = TRUE)
max(FullT$Pool.7, na.rm = TRUE)
max(FullT$Pool.8, na.rm = TRUE) #32.536
max(FullT$Pool.9, na.rm = TRUE)
max(FullT$Pool.10, na.rm = TRUE)
max(FullT$Pool.11, na.rm = TRUE)
max(FullT$Pool.13, na.rm = TRUE)
max(FullT$Pool.15, na.rm = TRUE) 
max(FullT$Pool.16, na.rm = TRUE)
max(FullT$Pool.27, na.rm = TRUE)
max(FullT$Pool.31, na.rm = TRUE)
max(FullT$Pool.33, na.rm = TRUE)
max(FullT$Pool.35, na.rm = TRUE) 
max(FullT$Pool.36, na.rm = TRUE)

# REMOVING NAs - NOT PRETTY but not in the mood to fight R-------------------------------------------------
FullT <-FullT %>% filter(!is.na(Pool.5))
FullT <-FullT %>% filter(!is.na(Pool.7))
FullT <-FullT %>% filter(!is.na(Pool.8))
FullT <-FullT %>% filter(!is.na(Pool.9))
FullT <-FullT %>% filter(!is.na(Pool.10))
FullT <-FullT %>% filter(!is.na(Pool.11))
FullT <-FullT %>% filter(!is.na(Pool.13))
FullT <-FullT %>% filter(!is.na(Pool.15))
FullT <-FullT %>% filter(!is.na(Pool.16))
FullT <-FullT %>% filter(!is.na(Pool.27))
FullT <-FullT %>% filter(!is.na(Pool.31))
FullT <-FullT %>% filter(!is.na(Pool.33))
FullT <-FullT %>% filter(!is.na(Pool.35))
FullT <-FullT %>% filter(!is.na(Pool.36))

DayTotalT<-FullT %>% filter(DayNight!="Night")
NightTotalT<-FullT %>% filter(DayNight!="Day")


#Remove any NAs (should be none)------
#DayTotalT <- na.omit(DayTotalT)


#Deal with formatting
#Only need below for detided temps csv
#DayTotalT$Date <- sub("^0+", "", DayTotalT$Date)
# Now, convert the Date column to a Date object
#Pay attention is data is - or / between dates. Change as needed
DayTotalT$Date <- as.Date(DayTotalT$Date, format="%m/%d/%y")

#Break Temp data up by Treatment so you can see better-----

#CO2 + Warmed------------
par(mgp = c(2.6, 1, 0))
plot(DayTotalT$Pool.5~Date, lwd=3, cex.axis=1.5, cex.lab=1.5, ylim=c(-1.5,33), xlab="Date", ylab="Temperature (°C)", type="l", lty=3, col = rgb(178/255, 34/255, 34/255, alpha = 0.3), data=DayTotalT)

lines(DayTotalT$Pool.11~Date, type="l", lty=3, lwd=3, col = rgb(178/255, 34/255, 34/255, alpha = 0.3), data=DayTotalT)
lines(DayTotalT$Pool.36~Date, type="l",  lty=3, lwd=3, col = rgb(178/255, 34/255, 34/255, alpha = 0.3), data=DayTotalT)

specific_date <- as.Date("2019-04-01")  # Specify the date
abline(v = specific_date, 
       col = rgb(0, 0, 0, alpha = 0.35),  #with 35% opacity
       lwd = 2, 
       lty = 2) 

#Warmed------------
par(mgp = c(2.6, 1, 0))
plot(DayTotalT$Pool.8~Date, lwd=3, cex.axis=1.5, cex.lab=1.5, ylim=c(-1.5,33), xlab="Date", ylab="Temperature (°C)", type="l", lty=1, col = rgb(178/255, 34/255, 34/255, alpha = 0.3), data=DayTotalT)

lines(DayTotalT$Pool.13~Date, type="l", lty=1, lwd=3, col = rgb(178/255, 34/255, 34/255, alpha = 0.3), data=DayTotalT)
lines(DayTotalT$Pool.16~Date, type="l", lty=1, lwd=3, col = rgb(178/255, 34/255, 34/255, alpha = 0.3), data=DayTotalT)
lines(DayTotalT$Pool.33~Date, type="l", lty=1, lwd=3, col = rgb(178/255, 34/255, 34/255, alpha = 0.3), data=DayTotalT)

specific_date <- as.Date("2019-04-01")  # Specify the date
abline(v = specific_date, 
       col = rgb(0, 0, 0, alpha = 0.35),  #with 35% opacity
       lwd = 2, 
       lty = 2) 

#Unmanipulated------------
par(mgp = c(2.6, 1, 0))
plot(DayTotalT$Pool.7~Date, lwd=3, cex.axis=1.5, cex.lab=1.5, ylim=c(-1.5,33), xlab="Date", ylab="Temperature (°C)", type="l", lty=1, col = rgb(120/255, 120/255, 120/255, alpha = 0.3), data=DayTotalT)

lines(DayTotalT$Pool.10~Date, type="l", lty=1, lwd=3, col = rgb(120/255, 120/255, 120/255, alpha = 0.3), data=DayTotalT)
lines(DayTotalT$Pool.31~Date, type="l", lty=1, lwd=3, col = rgb(120/255, 120/255, 120/255, alpha = 0.3), data=DayTotalT)

specific_date <- as.Date("2019-04-01")  # Specify the date
abline(v = specific_date, 
       col = rgb(0, 0, 0, alpha = 0.35),  #with 35% opacity
       lwd = 2, 
       lty = 2) 

#CO2 Added------------
par(mgp = c(2.6, 1, 0))
plot(DayTotalT$Pool.9~Date, lwd=3, cex.axis=1.5, cex.lab=1.5, ylim=c(-1.5,33), xlab="Date", ylab="Temperature (°C)", type="l", lty=3, col = rgb(120/255, 120/255, 120/255, alpha = 0.3), data=DayTotalT)

lines(DayTotalT$Pool.15~Date, type="l", lty=3, lwd=3, col = rgb(120/255, 120/255, 120/255, alpha = 0.3), data=DayTotalT)
lines(DayTotalT$Pool.27~Date, type="l", lty=3, lwd=3, col = rgb(120/255, 120/255, 120/255, alpha = 0.3), data=DayTotalT)
lines(DayTotalT$Pool.35~Date, type="l", lty=3, lwd=3, col = rgb(120/255, 120/255, 120/255, alpha = 0.3), data=DayTotalT)

specific_date <- as.Date("2019-04-01")  # Specify the date
abline(v = specific_date, 
       col = rgb(0, 0, 0, alpha = 0.35),  # with 35% opacity
       lwd = 2, 
       lty = 2) 

####NIGHTTIME TEMPERATURE FIGURE- Supplementary Figure X-------------------------------------------------------------------------
setwd("/Users/racinerangel/Desktop/Temperature")
#NightTotalT<-read.csv("FullNighttimeLowTempCL.csv")

#Remove any NAs (should be none)------
NightTotalT <- na.omit(NightTotalT)

#Deal with formatting
#NightTotalT$Date <- sub("^0+", "", NightTotalT$Date)
# Now, convert the Date column to a Date object
#NightTotalT$Date <- as.Date(NightTotalT$Date, format="%Y-%m-%d")

NightTotalT$Date <- as.Date(NightTotalT$Date, format="%m/%d/%y")

#Break Temp data up by Treatment so you can see better-----

#CO2 + Warmed------------
par(mgp = c(2.6, 1, 0))
plot(NightTotalT$Pool.5~Date, lwd=3, cex.axis=1.5, cex.lab=1.5, ylim=c(-1.5,33), xlab="Date", ylab="Temperature (°C)", type="l", lty=3, col = rgb(178/255, 34/255, 34/255, alpha = 0.3), data=NightTotalT)

lines(NightTotalT$Pool.11~Date, type="l", lty=3, lwd=3, col = rgb(178/255, 34/255, 34/255, alpha = 0.3), data=NightTotalT)
lines(NightTotalT$Pool.36~Date, type="l",  lty=3, lwd=3, col = rgb(178/255, 34/255, 34/255, alpha = 0.3), data=NightTotalT)

specific_date <- as.Date("2019-04-01")  # Specify the date
abline(v = specific_date, 
       col = rgb(0, 0, 0, alpha = 0.35),  #with 35% opacity
       lwd = 2, 
       lty = 2) 

#Warmed------------
par(mgp = c(2.6, 1, 0))
plot(NightTotalT$Pool.8~Date, lwd=3, cex.axis=1.5, cex.lab=1.5, ylim=c(-1.5,33), xlab="Date", ylab="Temperature (°C)", type="l", lty=1, col = rgb(178/255, 34/255, 34/255, alpha = 0.3), data=NightTotalT)

lines(NightTotalT$Pool.13~Date, type="l", lty=1, lwd=3, col = rgb(178/255, 34/255, 34/255, alpha = 0.3), data=NightTotalT)
lines(NightTotalT$Pool.16~Date, type="l", lty=1, lwd=3, col = rgb(178/255, 34/255, 34/255, alpha = 0.3), data=NightTotalT)
lines(NightTotalT$Pool.33~Date, type="l", lty=1, lwd=3, col = rgb(178/255, 34/255, 34/255, alpha = 0.3), data=NightTotalT)

specific_date <- as.Date("2019-04-01")  # Specify the date
abline(v = specific_date, 
       col = rgb(0, 0, 0, alpha = 0.35),  #with 35% opacity
       lwd = 2, 
       lty = 2) 

#Unmanipulated------------
par(mgp = c(2.6, 1, 0))
plot(NightTotalT$Pool.7~Date, lwd=3, cex.axis=1.5, cex.lab=1.5, ylim=c(-1.5,33), xlab="Date", ylab="Temperature (°C)", type="l", lty=1, col = rgb(120/255, 120/255, 120/255, alpha = 0.3), data=NightTotalT)

lines(NightTotalT$Pool.10~Date, type="l", lty=1, lwd=3, col = rgb(120/255, 120/255, 120/255, alpha = 0.3), data=NightTotalT)
lines(NightTotalT$Pool.31~Date, type="l", lty=1, lwd=3, col = rgb(120/255, 120/255, 120/255, alpha = 0.3), data=NightTotalT)

specific_date <- as.Date("2019-04-01")  # Specify the date
abline(v = specific_date, 
       col = rgb(0, 0, 0, alpha = 0.35),  # with 35% opacity
       lwd = 2, 
       lty = 2) 

#CO2 Added------------
par(mgp = c(2.6, 1, 0))
plot(NightTotalT$Pool.9~Date, lwd=3, cex.axis=1.5, cex.lab=1.5, ylim=c(-1.5,33), xlab="Date", ylab="Temperature (°C)", type="l", lty=3, col = rgb(120/255, 120/255, 120/255, alpha = 0.3), data=NightTotalT)

lines(NightTotalT$Pool.15~Date, type="l", lty=3, lwd=3, col = rgb(120/255, 120/255, 120/255, alpha = 0.3), data=NightTotalT)
lines(NightTotalT$Pool.27~Date, type="l", lty=3, lwd=3, col = rgb(120/255, 120/255, 120/255, alpha = 0.3), data=NightTotalT)
lines(NightTotalT$Pool.35~Date, type="l", lty=3, lwd=3, col = rgb(120/255, 120/255, 120/255, alpha = 0.3), data=NightTotalT)

specific_date <- as.Date("2019-04-01")  # Specify the date
abline(v = specific_date, 
       col = rgb(0, 0, 0, alpha = 0.35),  #with 35% opacity
       lwd = 2, 
       lty = 2) 