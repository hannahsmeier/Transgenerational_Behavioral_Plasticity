#Code for figure 2, 4 and S6 in Meier et al., 2022

library(tidyverse)
library(gridExtra)
library(mgcv) #GAM
install.packages("bbmle")
library(bbmle) #has AICtab function
library(car)
library(dplyr)
library(ggplot2)

#######################
#Initial Set Up#######
#######################


#load data
#paths for data have been removed but data file names remain
dat <- read_csv("results/combined_bootstrap_results_v2.csv")

speed_density <- read_csv("~/Desktop/thesis_data/speed_density.csv") %>%
  select(date,treatment, day,RFU)

#Dat2 includes RFU
dat2 <- merge(dat, speed_density, by=c("treatment", "day")) 

#################################################
#Day 0 Bar Graphs used to create figure 2#
#################################################

dat_sum2<-dat2 %>% 
  filter(day=="0")%>%
  group_by (treatment,Acc_Temp,Acute_Temp, day) %>%
  summarise(per_mov=mean(pmove),
            speed=mean(med),
            angle=mean(mu),
            mean_RFU=mean(RFU)) 

#Filtering only for day 0
#Averaging across treatments
day0_dat <- dat_sum2 

day0_dat1 <- day0_dat %>%
  group_by (Acc_Temp, day) %>%
  summarise(mov=mean(per_mov),
            speed_mean=mean(speed),
            angle_mean=mean(angle))
day0 <- dat %>%
  filter(day=="0")

#Running one-way ANOVAs each response variable only on day 0
#This is using the raw data looking at acclimation temperature effects

#Here testing whether RFU has a significant effect on any variables
rmaovmod_FE_0_speed<-aov(data= dat_sum2 %>% filter (day==0) , speed ~ as.factor(Acc_Temp)+mean_RFU)
summary(rmaovmod_FE_0_speed)

rmaovmod_FE_0_mov<-aov(data= dat_sum2 %>% filter (day==0) , per_mov ~  as.factor(Acc_Temp)+mean_RFU)
summary(rmaovmod_FE_0_mov)

rmaovmod_FE_0_angle<-aov(data= dat_sum2 %>% filter (day==0) , angle ~ as.factor(Acc_Temp)+mean_RFU)
summary(rmaovmod_FE_0_angle)

#Here to see whether acclimation temperture is significant on day 0
rmaovmod_FE_0_speed<-aov(data= dat %>% filter (day==0) , med ~ as.factor(Acc_Temp))
summary(rmaovmod_FE_0_speed)

rmaovmod_FE_0_mov<-aov(data= dat %>% filter (day==0) , pmove ~  as.factor(Acc_Temp))
summary(rmaovmod_FE_0_mov)

rmaovmod_FE_0_angle<-aov(data= dat  %>% filter (day==0) , mu ~  as.factor(Acc_Temp))
summary(rmaovmod_FE_0_angle)

#Using Holm-bonferroni on day 0 acclimation temperature significance

day0_p <- c(1.81e-12, 0.0576, 6.38e-05)

p<- p.adjust(day0_p, "holm")
print(p)

#Summary stats for day 0
day0_sum <- dat %>%
  filter(day == "0") %>%
  filter(Acc_Temp == Acute_Temp) %>%
  group_by(Acc_Temp, Acute_Temp) %>%
  summarise(mean_speed = mean(med),
            sd_speed=sd(med),
            mean_move = mean(pmove),
            sd_move=sd(pmove),
            mean_angle = mean(mu),
            sd_angle= sd(mu)) %>%
  mutate(speed_us = mean_speed*20) %>%
  mutate(speed_us_sd = sd_speed*20) %>%
  mutate(speed_se = speed_us_sd/sqrt(length((speed_us))),
         angle_se = sd_angle/sqrt(length(mean_angle)),
         move_se = sd_move/sqrt(length(mean_move)))

#Building a bar plot for the effects of acclimation temperature on each response variable
day0_speed_bar<-ggplot(day0_sum, mapping =  aes(fill=as.factor(Acc_Temp),x=as.factor(Acc_Temp), y=speed_us)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=speed_us-speed_se,ymax=speed_us+speed_se),width=.2, position=position_dodge(.9))+
  scale_fill_manual(values = c("blue", "orange","red"))+
  labs(fill="Acclimated Temp (°C)")+
  scale_y_continuous("Median Speed")+
  scale_x_discrete("Acclimated Temp (°C)",breaks=c(12.5, 25, 37.5),labels=c(12.5, 25, 37.5)) +
  theme_bw(base_size=18) +
  theme(legend.position = "None") +
  theme(aspect.ratio = 1)


day0_speed_bar

day0_mov_bar<-ggplot(day0_sum, mapping =  aes(fill=as.factor(Acc_Temp),x=as.factor(Acc_Temp), y=mean_move)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean_move-move_se,ymax=mean_move+move_se), width=.2, position=position_dodge(.9))+
  scale_fill_manual(values = c("blue", "orange","red"))+
  labs(fill="Acclimated Temp (°C)")+
  scale_y_continuous("Proportion Moving")+
  scale_x_discrete("Acclimated Temp (°C)",breaks=c(12.5, 25, 37.5),labels=c(12.5, 25, 37.5)) +
  theme_bw(base_size=18) +
  theme(legend.position="None") +
  #theme(strip.background = element_blank(),
  #      strip.text.x = element_blank(),
  #      legend.position = "None") +
  theme(aspect.ratio = 1)

day0_mov_bar

day0_angle_bar<-ggplot(day0_sum, mapping =  aes(fill=as.factor(Acc_Temp),x=as.factor(Acc_Temp), y=mean_angle)) +
  geom_bar(stat="identity")  +
  geom_errorbar(aes(ymin=mean_angle-angle_se,ymax=mean_angle+angle_se), width=.2, position=position_dodge(.9))+
  scale_fill_manual(values = c("blue", "orange","red"))+
  labs(fill="Acclimated Temp (°C)")+
  scale_x_discrete("Acclimated Temp (°C)",breaks=c(12.5, 25, 37.5),labels=c(12.5, 25, 37.5)) +
  scale_y_continuous("Average Absolute Angle") +
  theme_bw(base_size=18) +
  theme(legend.position="None") +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())  +
  theme(aspect.ratio = 1)

day0_angle_bar

#Combine these graphs
grid.arrange(day0_mov_bar,day0_speed_bar,day0_angle_bar,nrow =1)



#################################################
#CTA for figure 4 a and b#
#################################################

require(devtools)
#install_version("vegclust", version = "1.7.7", repos = "http://cran.us.r-project.org")
library(vegclust)
library(vegan)
library(RColorBrewer)
library(smacof)
#install.packages("ecotraj")
library(ecotraj)


#examine distributions of each variable
dat1 <- dat %>%
  mutate(speed_us = med*20)  %>%
  group_by(Acc_Temp, Acute_Temp, day) %>%
  summarise(speed = mean(speed_us),
            angle=mean(mu),
            per_mov = mean(pmove))
hist(dat1$angle) 
hist(log(dat1$speed)) 
hist(dat1$per_mov) 


#######
#####Now look at convergence plots
######

par(mar=c(5,5.5,1,1),oma=c(.5,.5,.5,0),mfrow=c(3,1),cex.axis=3,cex.lab=3)

colnames(dat1)
tmat<- dat1 %>% select("Acc_Temp", "Acute_Temp", "day", 
                          "per_mov", "speed", "angle") %>%
  filter(Acute_Temp==12.5) %>%
  mutate(logspeed=log(speed))


#Creating matrix of traits across time
C12.5= dist(tmat[c(4,6,7)])
#Distance matrix
C12.5= vegdist(tmat[c(4,6,7)],method="bray")

#Arrows represent species scores 
#Axis tells us % of variance explained by each axis 
#So acclimation temp. is more important for explaining patterns than day of trial
trajectoryPCoA(C12.5, tmat$Acc_Temp, tmat$day, lwd=0,plot=FALSE,traj.colors = "white")
grid(nx = 4, ny = 4, col = "lightgray", lwd = 2, lty=1, equilogs = TRUE)
grid(nx = 8, ny = 8, col = "lightgray", lwd = .25, lty=1, equilogs = TRUE)
par(new=TRUE)
trajectoryPCoA(C12.5, tmat$Acc_Temp, tmat$day, traj.colors = c("Blue","Orange", "Red"), lwd = 2,  survey.labels = F)
#legend("topright", legend=expression("Acute = 12.5"^O*" C"), bty="n", cex=3)
trajectoryConvergence(C12.5, tmat$Acc_Temp, tmat$day, symmetric = TRUE)

sp.C1<-trajectoryLengths(C12.5, tmat$Acc_Temp, tmat$day)


sp.C <-as.data.frame(t(sp.C1)) %>% dplyr::slice(1:(n()-3)) %>%  mutate(Acute_Temp=12.5) %>%
  pivot_longer(1:3,values_to = "Length", names_to="Acc_Temp") 


##now 25C
tmat<- dat1 %>% select("Acc_Temp", "Acute_Temp", "day", 
                          "per_mov", "speed", "angle") %>%
  filter(Acute_Temp==25) %>%
  mutate(logspeed=log(speed))

C25= dist(tmat[c(4,6,7)])
C25= vegdist(tmat[c(4,6,7)],method="bray")

#Arrows represent species scores 
#Axis tells us % of variance explained by each axis 
#So acclimation temp. is more important for explaining patterns than day of trial
trajectoryPCoA(C25, tmat$Acc_Temp, tmat$day, plot=FALSE,traj.colors = "white")
grid(nx = 4, ny = 4, col = "lightgray", lwd = 2, lty=1, equilogs = TRUE)
grid(nx = 8, ny = 8, col = "lightgray", lwd = .25, lty=1, equilogs = TRUE)
par(new=TRUE)
trajectoryPCoA(C25, tmat$Acc_Temp, tmat$day, traj.colors = c("Blue","Orange", "Red"), lwd = 2, xlim=c(-6,2),  survey.labels = F)
#legend("topright", legend=expression("Acute = 25"^O*" C"), bty="n", cex=3)
trajectoryConvergence(C25, tmat$Acc_Temp, tmat$day, symmetric = TRUE)

sp.W1<-trajectoryLengths(C25, tmat$Acc_Temp, tmat$day)

sp.W <-as.data.frame(t(sp.W1)) %>% dplyr::slice(1:(n()-3)) %>%  mutate(Acute_Temp=25) %>%
  pivot_longer(1:3,values_to = "Length", names_to="Acc_Temp") 

#37.5C
tmat<- dat1 %>% select("Acc_Temp", "Acute_Temp", "day",  
                          "per_mov", "speed", "angle") %>%
  filter(Acute_Temp==37.5) %>%
  mutate(logspeed=log(speed)) %>%
  drop_na()

C37.5= dist(tmat[c(4,6,7)])
C37.5= vegdist(tmat[c(4,6,7)],method="bray")

#Arrows represent species scores 
#Axis tells us % of variance explained by each axis 
#So acclimation temp. is more important for explaining patterns than day of trial
trajectoryPCoA(C37.5, tmat$Acc_Temp, tmat$day, plot=FALSE,traj.colors = "white")
grid(nx = 4, ny = 4, col = "lightgray", lwd = 2, lty=1, equilogs = TRUE)
grid(nx = 8, ny = 8, col = "lightgray", lwd = .25, lty=1, equilogs = TRUE)
par(new=TRUE)
trajectoryPCoA(C37.5, tmat$Acc_Temp, tmat$day, traj.colors = c("Blue","Orange", "Red"), lwd = 2,  survey.labels = F)

#legend("topright", col="White",  legend=expression("Acute = 37.5"^O*" C"), bty="n", lty=1, lwd = 0,cex=3)


trajectoryConvergence(C37.5, tmat$Acc_Temp, tmat$day, symmetric = TRUE)
sp.H1<-trajectoryLengths(C37.5, tmat$Acc_Temp, tmat$day)


sp.H <-as.data.frame(t(sp.H1)) %>% dplyr::slice(1:(n()-3)) %>%  mutate(Acute_Temp=37.5) %>%
  pivot_longer(1:3,values_to = "Length", names_to="Acc_Temp") 

#combine path lengths
sp.T<-rbind(sp.C,sp.W,sp.H)
sp.T$Acc_Temp<-as.factor(sp.T$Acc_Temp)
sp.T$Acute_Temp<-as.factor(sp.T$Acute_Temp)

sp.T_Sum<-sp.T%>%group_by(Acute_Temp,Acc_Temp)%>%summarise(Mean_Length=mean(Length),SE_Length=sd(Length)/sqrt(nrow(sp.T)/9))

#plot path lengths
path_length <- ggplot(data=sp.T_Sum,aes(x=Acute_Temp,y=Mean_Length,fill=Acc_Temp))+
  geom_bar(stat="identity",color= "black", position= position_dodge(.9))+ 
  geom_errorbar(aes(ymin=Mean_Length-SE_Length,ymax=Mean_Length+SE_Length), width=.2, position=position_dodge(.9))+
  scale_x_discrete("Acute Temperature")+
  scale_y_continuous("Mean Pathlength")+
  theme(legend.position="top") +
  scale_fill_manual(values=c("Blue","Orange","Red")) +
  theme_bw(base_size=18)
path_length

grid.arrange(cold, am, hot, path_length, nrow=2)

#test
Pmod<-lm(Length~ Acute_Temp * Acc_Temp, dat1a= sp.T )
anova(Pmod)

#export as pdf fig

#################################################
#Investigating the relationship between speed 
#and density for figure S6#
#################################################


##Add scaled speed
speed_density <- dat2 %>%
  mutate(speed_us= med*20) 

#Create summary statistic
speed_density_mean <- speed_density %>%
  group_by(Acc_Temp, Acute_Temp, day) %>%
  summarise(mean_RFU = mean(RFU),
            speed = mean(speed_us))



install.packages("ggpmisc")
library(ggpmisc)
library(tidyverse)
library(gridExtra)
library(mgcv) #GAM
library(bbmle) #has AICtab function
library(car)
library(dplyr)
library(ggplot2)

library(tidyverse)
library(devtools)
#install.packages("ggpubr")
library(ggpubr)

#Correct RFU to log10 count
speed_density_mean1 <- speed_density_mean %>%
  mutate(Count = ((11.54*mean_RFU)-566.8)) %>%
  mutate(Count_log = log10(Count))

#This is just to aid with the facet labels
speed_density_mean2 <- speed_density_mean1 %>%
  mutate(acc_temp = ifelse(Acc_Temp == "12.5", "cold",
                           ifelse(Acc_Temp == "25", "ambient", "hot"))) 
#Added to get the facets in the right order
speed_density_mean2$acc_order = factor(speed_density_mean2$acc_temp, levels=c("cold", "ambient","hot"))

#Graphing this model and creating stats
ggscatter(speed_density_mean2, x = "Count_log", y = "speed", 
          facet.by = c("Acute_Temp","acc_order"), levels=c("cold", "ambient", "hot"), 
          color ="acc_order",
          palette = c(cold="Blue",ambient="Orange", hot="Red"),
          add = "reg.line", conf.int = TRUE, 
          xlab = "Log10 Density (Cells/mL)", ylab = "Median Speed (microns/sec)") +
  theme_bw(base_size=18) +
  stat_cor(aes(color = acc_temp, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), show.legend = TRUE, size =8) +
  scale_fill_manual(values = c("Blue","Orange", "Red"),name = "Acclimated Temp (°C)",labels = c("cold", "ambient","hot") 
  ) +
  scale_colour_manual(values = c("Blue","Orange", "Red"),name = "Acclimated Temp (°C)", labels = c("cold", "ambient","hot")) 


