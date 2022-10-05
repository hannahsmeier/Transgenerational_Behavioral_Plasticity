
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Code for analyzing movement behavior of C. reinhardtii

# Code used for creating figure 3 Meier, et al. 2022

# CTK, September 2022

# Formerly "gamma_bootstrapping_v6.R"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#### Outline of analyses: ####

# (0) Preliminaries
# (1) Illustrate time series of individual cells
# (2) Movement analysis
#     A. Characterize velocity of heat-killed cells
#     B. Determine when cells are moving and characterize median velocity 
#        and how often cells move (point estimates)
#     C. Bootstrap estimate of confidence interval (median velocity)
#     D. Bootstrap estimate of confidence interval (proportion moving)
#     E. Point estimate & bootstrapping (turning angle)
#     F. Combine all point estimates and confidence intervals
# (3) Analyze differences in movement metrics across treamtents
#     A. Velocity
#     B. Proportion moving
#     C. Turning angle

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


#### (0) Preliminaries ####

# Load libraries:
library(ggplot2)
library(dplyr)
library(scales)
library(gamlss)
library(mgcv)
library(bbmle)
library(gridExtra)
library(bde)
library(mleTools)
#install.packages("betareg")
library(betareg)

# Define supplemental functions:
checkMove<-function(x) as.numeric(max(x)>0)
expit<-function(x){exp(x)/(1+exp(x))}
logit<-function(x){log(x/(1-x))}	

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#### (1) Illustrate time series of individual cells ####

# Load data:
idat<-read.csv("/Users/colin/Research/Active/Acclimation/movement/data/Links in tracks statistics.csv")
head(idat)

# For some exploratory plotting:
uids<-unique(idat$TRACK_ID)
plotTS<-function(i){
  tmp<-idat[idat$TRACK_ID==uids[i],]
  plot(tmp$VELOCITY,type='l',xlab='Time interval',ylab='Velocity')
}

# Pull out instance number 26:
plotTS(26)

tmp<-idat[idat$TRACK_ID==uids[26],]
len<-nrow(tmp)
kval<-min(c(floor(5+((80-5)/800)*len),80))

# GAM approach:
gam1<-gam(VELOCITY~s(EDGE_TIME,k=kval),
          data=tmp,family=Gamma(link="log"))

# generate vector of predictions based on this fit
pds<-predict(gam1,type='response',se.fit=T)
tmp$V.fit<-pds$fit
tmp$V.lwr95<-pds$fit-1.96*pds$se.fit  # calculate confidence band

head(tmp)

p1<-ggplot(tmp,aes(x=EDGE_TIME,y=VELOCITY))+
  geom_line(color='gray')+
  geom_line(aes(y=V.fit),color='black')+
  geom_line(aes(y=V.lwr95),color='black',linetype=2)+
  geom_hline(yintercept=5*0.32029244,color='red')+
  theme_bw()

p2<-tmp %>% 
  ggplot()+
  geom_histogram(aes(x=VELOCITY))+
  coord_cartesian(xlim=c(0,10))+
  theme_bw()+
  ggtitle("Full distribution")
p2

p3<-tmp %>% filter(V.lwr95>5*0.32029244) %>%
  ggplot()+
  geom_histogram(aes(x=VELOCITY))+
  coord_cartesian(xlim=c(0,10))+
  theme_bw()+
  ggtitle("Velocities while consistently moving")
p3

lay1<-rbind(c(1,1,1,1),
            c(2,2,3,3))
grid.arrange(p1,p2,p3,nrow=2,layout_matrix=lay1)

a1<-arrangeGrob(p1,p2,p3,nrow=2,layout_matrix=lay1)

ggsave("./results/illustrative_time_series.pdf",a1,width=11,height=7)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#### (2) Movement analysis ####

####    A. Characterize velocity of heat-killed cells ####

# load data:
ndat<-read.csv("./data/null hypotheses updated/Track statistics.csv")
head(ndat)

# truncate minimum velocity to avoid errors caused by exact zeros:
ndat$TRACK_MEAN_SPEED2<-ifelse(ndat$TRACK_MEAN_SPEED==0,0.01,ndat$TRACK_MEAN_SPEED)

# this influences less than 1% of the data...
1-length(ndat$TRACK_MEAN_SPEED[ndat$TRACK_MEAN_SPEED>0])/nrow(ndat)

# Fit gamma distribution to variation in velocity for these heat-killed particles
nd1<-mle2(TRACK_MEAN_SPEED2~dgamma2(mu=mu1,s=s),start=list(mu1=(0.3),s=0.03),data=ndat)
summary(nd1)

# compute likelihood profiles given model fit:
pf<-profile(nd1)

# Use likelihood profile to compute 99% confidence interval for the mean of the gamma distribution:
confint(pf,level = 0.99)
#     0.5 %     99.5 %
# mu1 0.27073745 0.32029244
# s   0.02008659 0.03421055

# This value is the basis for the threshold we use to define minimum velocity for a meaningfully moving cell,
#5*0.32029244


####    B. Determine when cells are moving and characterize median velocity ####
#        and how often cells move (point estimates)

#### Data import process ####

# Fix the root directory
root<-"/Users/colin/Research/Active/Acclimation/movement/"
setwd(root)

# grab high level fnames
setwd("./data/phyto_mov_data/")

# this creates a list of all files (with their file paths) in all subdirectories and saves it as a csv
system("du -a . > all_fnames2.csv")

# restore working directory
setwd(root)

# absorb the csv
ack<-read.csv("./data/phyto_mov_data/all_fnames2.csv",header = F)

# thin files to look only at the "links" files
ack2<-ack[grepl("links_in_tracks_statistics.csv",ack[,1]),]

# trim off some junk
ack3<-sapply(ack2,FUN=function(x) strsplit(x,split="\t")[[1]][2])
attr(ack3,"names")<-NULL

# drop 'old' folders
ack3<-ack3[!grepl("old",ack3)]

# extract dates and treatments from file names
ack4<-sapply(ack3,FUN=function(x) strsplit(x,split="/")[[1]][4])
dates<-sapply(ack4,FUN=function(x) strsplit(x,split="_")[[1]][1])
attr(dates,"names")<-NULL
treatments<-sapply(ack4,FUN=function(x) strsplit(x,split="_")[[1]][2])
attr(treatments,"names")<-NULL

# combine to create metadata table:
metadata<-data.frame(fpaths=ack3,dates,treatments)
head(metadata)

# load auxiliary data for creating date table:
datF<-read.csv("./data/FIJI_fin_summary.csv")
lookup<-unique(datF[,c('date','day')])
names(lookup)[1]<-"dates"

# augment metadata table with numeric 'day' code:
metadata<-merge(metadata,lookup,all.x=T)

# add more complete info on path
metadata$fpaths<-gsub("^.","./data/phyto_mov_data",metadata$fpaths)
#head(metadata)

# Now we can cycle through all of the links files, and conduct an analysis to determine which cells are moving, and when.


#### Determine movement status ####

# Define minimum speed used to declare a cell is moving
threshold<-5*0.32029244

# Create storage for results
res1<-data.frame(treatment=rep(NA,nrow(metadata)),day=NA,mu=NA,sig=NA,med=NA,pmove=NA,cmove=NA)

# Cycle through all data sets:
for(i in 1:nrow(metadata)){
  print(i/nrow(metadata))
  
  # grab the name of data file
  data.set<-metadata[i,]
  
  # load data for this treatment/replicate
  idat<-read.csv(data.set$fpaths)
  uids<-unique(idat$TRACK_ID)
  
  ### Plot time series
  print("Initial plotting")
  
  plotname<-paste("tracksTS",data.set$dates,data.set$treatments,data.set$day,sep="_")

  # (only create time series plot if it doesn't already exist) 
  # turn this OFF if there is a change in methods/plotting options.
  if(!file.exists(paste("./results/trackTS/",plotname,".png",sep=""))){
    png(filename = paste("./results/trackTS/",plotname,".png",sep=""),width=2*480,height=480)
    plot(VELOCITY~EDGE_TIME,data=idat[idat$TRACK_ID==uids[1],],type='l',
         col=alpha("black",0.1),ylim=c(0,14),xlim=c(0,1200),main=paste(data.set$dates,data.set$treatments,data.set$day))
    for(j in 2:length(uids)){
      lines(VELOCITY~EDGE_TIME,data=idat[idat$TRACK_ID==uids[j],],col=alpha("black",0.1))
    }
    abline(0.29415958,0,col='blue')
    abline(5*0.32029244,0,col='red')
    dev.off()
  }else{
    print("TS plot skipped")
  }

  #### Filter velocities to obtain target data set:
  
  print("Filtering velocities")
  
  uids<-unique(idat$TRACK_ID)
  # drop zeros (incompatible with Gamma family GAM)
  idat$VELOCITY<-ifelse(idat$VELOCITY<0.0008936,0.0001,idat$VELOCITY)
  
  # set up new columns
  idat$VELOCITY2<-NA
  idat$moveQ<-NA
  
  #j<-sample(seq(1,length(uids),1),size = 1)
  for(j in 1:length(uids)){
    #print(j/length(uids))
    tmp<-idat[idat$TRACK_ID==uids[j],]
    
    len<-nrow(tmp)
    
    if(len>=10){
      
      kval<-min(c(floor(5+((80-5)/800)*len),80))
      
      # GAM approach:
      gam1<-gam(VELOCITY~s(EDGE_TIME,k=kval),
                data=tmp,family=Gamma(link="log"))
      
      # generate vector of predictions based on this fit
      pds<-predict(gam1,type='response',se.fit=T)
      tmp$V.fit<-pds$fit
      tmp$V.lwr95<-pds$fit-1.96*pds$se.fit  # calculate confidence band
      
      #decide what to keep:
      tmp$VELOCITY2<-ifelse(tmp$V.lwr95>5*0.32029244,
                            tmp$VELOCITY,NA)
      
      # fill these results back into the idat frame:
      idat$VELOCITY2[idat$TRACK_ID==uids[j]]<-	tmp$VELOCITY2
      idat$moveQ[idat$TRACK_ID==uids[j]]<-ifelse(tmp$V.lwr95>5*0.32029244,1,0)
    }# else = entries stay NA
  }
  
  # subset data to focus on cells that move:
  wdat<-idat %>% filter(moveQ==1)
  
  # check for no data case:
  if(nrow(wdat)==0){
    print("no moving cells!")
    nc<-ncol(wdat)
    wdat[i,]<-rep(NA,nc)
    wdat<-cbind(wdat,data.frame(treatment=NA,day=NA))
    wdat<-na.omit(wdat)
    
    # summary statistics (null)
    res1[i,]<-c(treatment=data.set$treatments,day=data.set$day,
                mu=NA,sig=NA,med=NA,pmove=0,cmove=0)
  }else{
    wdat$treatment<-data.set$treatments
    wdat$day<-data.set$day
    
    ### Calculations
    print("Characterizing distribution")
    
    # proportion of all observation intervals where a cell is moving
    pmove<-sum(idat$moveQ)/nrow(idat)
    
    # proportion of cells that move at least once
    cellsmove<-idat %>% group_by(TRACK_ID) %>% summarise(cMove=checkMove(moveQ)) %>% ungroup() %>% summarise(sum(cMove))
    cellsmove<-cellsmove/length(unique(idat$TRACK_ID))
    
    # Characterize original data on the distribution of movement velocities:
    
    # calculate the median velocition 
    md<-median(wdat$VELOCITY2)
    
    # fit a normal distribution
    f1<-gamlss(VELOCITY2~1,family=NO,data=wdat,control=gamlss.control(trace=F)) 
    
    # fit a Gamma distribution
    f2<-gamlss(VELOCITY2~1,family=GA,data=wdat,control=gamlss.control(trace=F))
    
    orig<-c(exp(f2$mu.coefficients),exp(f2$sigma.coefficients),md)
    
    # save summary statistics (e.g., median, movement fractions, and Gamma fit results)
    res1[i,]<-c(treatment=data.set$treatments,day=data.set$day,
                mu=orig[1],sig=orig[2],med=orig[3],pmove=pmove,cmove=cellsmove)
    
    # plot distribution fits
    plotname2<-paste("dists",data.set$dates,data.set$treatments,data.set$day,sep="_")
    png(filename = paste("./results/dists/",plotname2,".png",sep=""),width=480,height=480)
    dns<-density(wdat$VELOCITY2,bw = 0.2)
    plot(dns,main=paste(data.set$dates,data.set$treatments,data.set$day))
    curve(dGA(x,exp(f2$mu.coefficients),exp(f2$sigma.coefficients)),0,12,col='red',add=T)
    curve(dNO(x,f1$mu.coefficients,sigma=exp(f1$sigma.coefficients)),0,12,col='blue',add=T)
    segments(md,0,md,1.2*max(dns$y),lty=3)
    dev.off()
  }
  
  # save file w/all observation intervals and whether or not cells were moving:
  c1<-strsplit(data.set$fpaths,split='/')[[1]][5]
  write.csv(wdat,paste("./data/filtered_phyto_mov_data/filtered_",c1,".csv",sep=""),row.names=F)
  
  # save file w/summary statistics of movement
  write.csv(res1,"./results/summary_stats_phyto_mov_v1.csv",row.names=F)
}



####  C. Bootstrap estimate of confidence interval (median velocity) ####

# Rationale: From the prior analysis, we have point estimates (e.g. medians) characterizing the velocity of moving cells in all of our different treatments and replicates. Before determining whether these values vary in interesting or meaningful ways among our treatments, we also must characterize how certain these estimates are. For this reason, we use a boot-strapping approach to estimate confidence intervals around our point estimates, and then later use the width of these confidence intervals to account for varying levels of uncertainty in the value of point estimates across replicates/treatments.

# number of bootstrapping replicates to run
nreps<-1000  

# data storage:
res2<-data.frame(treatment=rep(NA,nrow(metadata)),day=NA,med.lwr=NA,med.upr=NA)

for(i in 1:nrow(metadata)){
  print(i/nrow(metadata))
  
  # read data file
  data.set<-metadata[i,]
  c1<-strsplit(data.set$fpaths,split='/')[[1]][5]
  wdat<-read.csv(paste("./data/filtered_phyto_mov_data/filtered_",c1,".csv",sep=""))

  # check for no data case:
  if(nrow(wdat)==0){
    print("no moving cells!")
    
    # save empty results (can't estimate a confidence interval)
    res2[i,]<-c(treatment=data.set$treatments,day=data.set$day,
                med.lwr=NA,med.upr=NA)
  }else{
    ### Run bootstrapping:
    print("Bootstrapping")
    
    # function will resample velocity data from the original distribution,
    # with replacement, and compute the new median for the resulting distribution
    bootfunc<-function(k){
      tmp2<-data.frame(VELOCITY2=sample(wdat$VELOCITY2,nrow(wdat),replace = T))
      median(tmp2$VELOCITY2)
    }
    
    # call function nreps number of times, compile results
    meds<-unlist(lapply(seq(1,nreps,1),bootfunc))
    
    # figure out 2.5th and 97.5th quantile of the bootstrap distribution of medians
    med.ci<-quantile(meds,probs=c(0.025,0.975))
    
    # save results for this treatment/replicate
    res2[i,]<-c(treatment=data.set$treatments,day=data.set$day,
                med.lwr=med.ci[1],med.upr=med.ci[2])
  }  

  # Save results:
  write.csv(res2,"./results/bootstrap_velocities_v5.csv",row.names=F)
}

# format results
res2$med.lwr<-as.numeric(res2$med.lwr)
res2$med.upr<-as.numeric(res2$med.upr)


####  D. Bootstrap estimate of confidence interval (proportion moving) ####

# number of bootstrapping replicates
nreps<-1000  

# data storage:
res4<-data.frame(treatment=rep(NA,nrow(metadata)),day=NA,pmove.lwr=NA,pmove.upr=NA)

i<-1
for(i in 1:nrow(metadata)){
  print(i/nrow(metadata))
  
  # read file (only intervals where cells were deemed to be moving)
  data.set<-metadata[i,]
  c1<-strsplit(data.set$fpaths,split='/')[[1]][5]
  wdat<-read.csv(paste("./data/filtered_phyto_mov_data/filtered_",c1,".csv",sep=""))
  
  # check for no data case:
  if(nrow(wdat)==0){
    print("no moving cells!")
    
    # save empty results
    res4[i,]<-c(treatment=data.set$treatments,day=data.set$day,
                pmove.lwr=NA,pmove.upr=NA)
  }else{
    ### Run bootstrapping:
    print("Bootstrapping")
    
    # load original data set (both moving & not moving), so we know how often each occurred
    idat<-read.csv(data.set$fpaths)
    nlinks<-nrow(idat)

    # recapitulate vector of appropriate number of 0's and 1's:
    tmp1<-c(rep(1,nrow(wdat)),rep(0,I(nlinks-nrow(wdat))))
    
    # this function will resample the proportion moving data and compute a new proportion
    # moving given the new distribution
    bootfunc<-function(k){
      tmp2<-sample(tmp1,nlinks,replace = T)
      sum(tmp2)/nlinks
    }
    
    # call function nreps number of times, compile results
    pmoves<-unlist(lapply(seq(1,nreps,1),bootfunc))
    
    # estimate CI based on 2.5th and 97.5th quantile of the bootstrapped distribution
    pmove.ci<-quantile(pmoves,probs=c(0.025,0.975))
    
    # store results for this treatment/replicate
    res4[i,]<-c(treatment=data.set$treatments,day=data.set$day,
                pmove.lwr=pmove.ci[1],pmove.upr=pmove.ci[2])
  }  
  
  # Save results:
  write.csv(res4,"./results/bootstrap_pmove_v1.csv",row.names=F)
}

# format columns
res4$pmove.lwr<-as.numeric(res4$pmove.lwr)
res4$pmove.upr<-as.numeric(res4$pmove.upr)



####  E. Point estimate & bootstrapping (turning angle) ####

# Turning angle data arises from a separate analysis using the HMM package, which produced a set of files (one per treatment and replicate) containing estimates of turning angle for each observation interval for each cell. To analyze these in a similar fashion to the velocity and proportion moving data, we have to load each file, match it to the results of our prior analyses (indicating when each cell is 'moving'), characterize the resulting distribution of turning angles (generating point estimates), and then quantify uncertainty in these point estimates (bootstrapping)

# fix root directory
root<-"/Users/colin/Research/Active/Acclimation/movement/"
setwd(root)

# bore into the folder containing the angle data sets:
setwd("./data/angle data/")

# this creates a list of all files in this folder
system("ls > angle_fnames.csv")

# restore working directory
setwd(root)

# absorb the csv of file names into R
ack<-read.csv("./data/angle data/angle_fnames.csv",header = F)
ack<-ack[-1,] # trim
head(ack)

# set up data storage (point estimates and confidence intervals):
res5<-data.frame(treatment=rep(NA,length(ack)),day=NA,mu=NA,mu.lwr=NA,mu.upr=NA,phi=NA,phi.lwr=NA,phi.upr=NA)


# number of bootstrap replicates to run:
nreps<-1000

## CAUTION - this next step takes *days* to run ##

# loop through treatments/replicates, loading angle data and making calculations
for(i in 1:length(ack)){
  print(i/length(ack))
  
  # load data:
  a1<-read.csv(paste("./data/angle data/",ack[i],sep=""))
  
  # load matching data with info on which observation intervals cells were moving:
  f1<-read.csv(paste("./data/filtered_phyto_mov_data/filtered_phyto_mov_",a1$treatment[1],"_",a1$date[1],".csv",sep=""))
  
  # check for no data case:
  if(nrow(f1)==0){
    print("no moving cells!")
    
    # save null results
    res5[i,]<-data.frame(treatment=a1$treatment[1],day=a1$day[1],mu=NA,mu.lwr=NA,mu.upr=NA,phi=NA,phi.lwr=NA,phi.upr=NA)
    
  }else{
    
    # merge data
    f1.2<-merge(a1,f1[,c('TRACK_ID','Label','EDGE_TIME','moveQ')],all.x=T)
    
    # create scaled angle column (absolute turning angle divided by pi)
    f1.2$scaled.angle<-abs(f1.2$angle)/pi
    
    # truncate to avoid errors with beta distribution when exact 0's or 1's are encountered
    f1.2$scaled.angle<-ifelse(f1.2$scaled.angle<=I(1e-7),I(1e-7),f1.2$scaled.angle)
    f1.2$scaled.angle<-ifelse(f1.2$scaled.angle>=I(1-1e-7),I(1-1e-7),f1.2$scaled.angle)
    
    # create point estimates for the scaled turning angle of moving cells
    # in this case, by characterizing the distribution as a beta distribution
    br1<-betareg(scaled.angle~1,data=f1.2[!is.na(f1.2$moveQ),])
    cf<-coef(br1)
    
    # save this fit for later visual inspection:
    pdf(file = paste("./results/angle_dist_fits/",a1$treatment[1],"_",a1$date[1],".pdf",sep=""))
    plot(bde(f1.2$scaled.angle[!is.na(f1.2$moveQ)], estimator="boundarykernel",b=0.025),main=paste(a1$treatment[1],"_",a1$date[1]," ",round(expit(cf[1]),3),sep=""))
    curve(dbeta2(x,mu=expit(cf[1]),phi = cf[2]),0,1,col='red',add=T)
    dev.off()
    
    # Now, run bootstrapping:
    
    # this function will resample scaled angle data and re-fit a beta distribution
    # then return the resulting coefficients for the 'new' distribution
    bootfunc<-function(k){
      tmp2<-data.frame(scaled.angle=sample(f1.2$scaled.angle[!is.na(f1.2$moveQ)],nrow(f1.2),replace = T))
      
      br1<-betareg(scaled.angle~1,data=tmp2)
      cf<-coef(br1)
      return(cf)
    }
    
    print("Bootstrapping:")
    
    # call the bootfunc function nreps number of times, compile results
    betas<-unlist(lapply(seq(1,nreps,1),bootfunc))
    
    # parse the mean (mu) and scale (phi) parameter of the beta distribution arising from
    # the bootstrapped fits:
    mus<-expit(betas[seq(1,length(betas),2)])
    phis<-betas[seq(2,length(betas),2)]
    
    # estimate 95% confidence intervals
    mu.ci<-quantile(mus,probs=c(0.025,0.975))
    phi.ci<-quantile(phis,probs=c(0.025,0.975))
    
    # save point estimates and bootstrapped confidence intervals for this treatment/replicate
    res5[i,]<-data.frame(treatment=a1$treatment[1],day=a1$day[1],mu=expit(cf[1]),mu.lwr=mu.ci[1],mu.upr=mu.ci[2],phi=cf[2],phi.lwr=phi.ci[1],phi.upr=phi.ci[1])
  }
  
  # Save overall data file
  write.csv(res5,"./results/angle_bootstrap_results_v2.csv",row.names=F)
}
#head(res5)


####  F. Combine all point estimates and confidence intervals ####

# load output from above computations of point estimates and CI's:
res1<-read.csv("./results/summary_stats_phyto_mov_v1.csv")
res2<-read.csv("./results/bootstrap_velocities_v5.csv")
res4<-read.csv("./results/bootstrap_pmove_v1.csv")
res5<-read.csv("./results/angle_bootstrap_results_v2.csv")

# drop the 'mu' and 'sig' columns from res1, which are based on parametric fits to the distribution of velocities, as we've decided to only use the 'median' for subsequent analyses.
res1<-res1[,which(!(names(res1) %in% c('mu','sig')))]


# Merge:

dim(res1)

# point estimates and velocity CI's
res.all<-merge(res1,res2,all.x=T)
dim(res.all)

# pmoving CI's
res.all<-merge(res.all,res4,all.x=T)
dim(res.all)

# angle data
res.all<-merge(res.all,res5,all.x=T)
dim(res.all)


# compute width of CI's
res.all$med.ci.width<-res.all$med.upr-res.all$med.lwr
res.all$pmove.ci.width<-res.all$pmove.upr-res.all$pmove.lwr
res.all$mu.ci.width<-res.all$mu.upr-res.all$mu.lwr

# create look-up table matching treatment label to acclimated and acute temperatures:
lk.tab<-data.frame(treatment = c("1A", "1B", "1C", "2A", "2B", "2C", 
                   "3A", "3B", "3C", "4A", "4B", "4C", "5A", "5B", "5C", "6A", "6B", 
                   "6C", "7A", "7B", "7C", "8A", "8B", "8C", "9A", "9B", "9C"), 
     Acc_Temp = c(12.5, 12.5, 12.5, 12.5, 12.5, 12.5, 12.5, 12.5, 
                  12.5, 25, 25, 25, 25, 25, 25, 25, 25, 25, 37.5, 37.5, 37.5, 
                  37.5, 37.5, 37.5, 37.5, 37.5, 37.5), 
     Acute_Temp = c(12.5,12.5, 12.5, 25, 25, 25, 37.5, 37.5, 37.5, 12.5, 12.5, 12.5, 
                  25, 25, 25, 37.5, 37.5, 37.5, 12.5, 12.5, 12.5, 25, 25, 25, 
                  37.5, 37.5, 37.5))


res.all<-merge(res.all,lk.tab)

write.csv(res.all,"./results/combined_bootstrap_results_v2.csv",row.names=F)




#### (3) Analyze differences in movement metrics across treamtents ####

# Load point estimates & bootstrapped CI's:
res3<-read.csv("./results/combined_bootstrap_results_v2.csv")
head(res3)

####     A. Velocity ####


#### 12 C Acute:  ####
bob.12<-res3[res3$Acute_Temp==12.5,]
bob.12$Acc_Temp<-as.factor(bob.12$Acc_Temp)

weights <- 1/bob.12$med.ci.width
weights <- weights/mean(weights)

gm.12.4<-gam(med~s(day,by=Acc_Temp,k=4)+Acc_Temp,family=Gamma(link="log"),data=bob.12,weights = weights)
gm.12.5<-gam(med~s(day,by=Acc_Temp,k=5)+Acc_Temp,family=Gamma(link="log"),data=bob.12,weights = weights)
gm.12.6<-gam(med~s(day,by=Acc_Temp,k=6)+Acc_Temp,family=Gamma(link="log"),data=bob.12,weights = weights)
gm.12.7<-gam(med~s(day,by=Acc_Temp,k=7)+Acc_Temp,family=Gamma(link="log"),data=bob.12,weights = weights)
gm.12.8<-gam(med~s(day,by=Acc_Temp,k=8)+Acc_Temp,family=Gamma(link="log"),data=bob.12,weights = weights)

# favors k = 6-8, but on visual inspection, over-fit
AICctab(gm.12.4,gm.12.5,gm.12.6,gm.12.7,gm.12.8)

# use k=6
gm.12<-gam(med~s(day,by=Acc_Temp,k=6)+Acc_Temp,family=Gamma(link="log"),data=bob.12,weights = weights)
summary(gm.12)

gm.12.0<-gam(med~s(day,k=6)+Acc_Temp,family=Gamma(link="log"),data=bob.12,weights = weights)
summary(gm.12.0)

# so, still support for distinct smooths
AICctab(gm.12,gm.12.0)


# Family: Gamma 
# Link function: log 
# 
# Formula:
#   med ~ s(day, by = Acc_Temp, k = 6) + Acc_Temp
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   1.66809    0.02490  66.984  < 2e-16 ***
#   Acc_Temp25   -0.06839    0.03564  -1.919   0.0599 .  
# Acc_Temp37.5 -0.21456    0.03746  -5.727 3.74e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F  p-value    
# s(day):Acc_Temp12.5 1.000  1.000 0.152    0.698    
# s(day):Acc_Temp25   4.660  4.936 9.402 1.54e-06 ***
#   s(day):Acc_Temp37.5 4.843  4.986 7.729 2.33e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.626   Deviance explained = 66.6%
# GCV = 0.019805  Scale est. = 0.015854  n = 72

bob.12 %>% filter(day %in% c(0,7,14)) %>%
  ggplot(aes(y=med))+
  geom_boxplot(aes(fill=factor(Acc_Temp)))+
  facet_wrap(~day)+
  scale_fill_manual('Acclimated\nTemperature',values=c("blue","orange","red"))+
  theme_bw()

an.12.0<-glm(med~Acc_Temp,family=Gamma(link = "log"),weights = weights[bob.12$day==0],data=bob.12[bob.12$day==0,])
summary(an.12.0)

an.12.7<-glm(med~Acc_Temp,family=Gamma(link = "log"),weights = weights[bob.12$day==7],data=bob.12[bob.12$day==7,])
summary(an.12.7)

an.12.14<-glm(med~Acc_Temp,family=Gamma(link = "log"),weights = weights[bob.12$day==14],data=bob.12[bob.12$day==14,])
summary(an.12.14)




#### 25 C Acute: ####
bob.25<-res3[res3$Acute_Temp==25,]
bob.25$Acc_Temp<-factor(bob.25$Acc_Temp,levels=c('25','12.5','37.5'))

weights <- 1/bob.25$med.ci.width
weights <- weights/mean(weights)

gm.25.4<-gam(med~s(day,by=Acc_Temp,k=4)+Acc_Temp,family=Gamma(link="log"),data=bob.25,weights = weights)
gm.25.5<-gam(med~s(day,by=Acc_Temp,k=5)+Acc_Temp,family=Gamma(link="log"),data=bob.25,weights = weights)
gm.25.6<-gam(med~s(day,by=Acc_Temp,k=6)+Acc_Temp,family=Gamma(link="log"),data=bob.25,weights = weights)
gm.25.7<-gam(med~s(day,by=Acc_Temp,k=7)+Acc_Temp,family=Gamma(link="log"),data=bob.25,weights = weights)
gm.25.8<-gam(med~s(day,by=Acc_Temp,k=8)+Acc_Temp,family=Gamma(link="log"),data=bob.25,weights = weights)

# favors k = 7 and 8, but on visual inspection, over-fit
AICctab(gm.25.4,gm.25.5,gm.25.6,gm.25.7,gm.25.8)

# use k=6
gm.25<-gam(med~s(day,by=Acc_Temp,k=6)+Acc_Temp,family=Gamma(link="log"),data=bob.25,weights = weights)
summary(gm.25)

gm.25.0<-gam(med~s(day,k=6)+Acc_Temp,family=Gamma(link="log"),data=bob.25,weights = weights)
summary(gm.25.0)

# so, still support for distinct smooths
AICctab(gm.25,gm.25.0)

# Family: Gamma 
# Link function: log 
# 
# Formula:
#   med ~ s(day, by = Acc_Temp, k = 6) + Acc_Temp
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   1.89734    0.02423  78.312  < 2e-16 ***
#   Acc_Temp12.5 -0.02812    0.03609  -0.779    0.439    
# Acc_Temp37.5 -0.25753    0.03461  -7.441 6.09e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df      F  p-value    
# s(day):Acc_Temp25   4.653  4.934 14.305  < 2e-16 ***
#   s(day):Acc_Temp12.5 3.393  3.991 11.123 1.42e-06 ***
#   s(day):Acc_Temp37.5 4.363  4.779  4.639   0.0017 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.751   Deviance explained = 78.5%
# GCV = 0.019423  Scale est. = 0.01433   n = 72


bob.25 %>% filter(day %in% c(0,7,14)) %>%
  ggplot(aes(y=med))+
  geom_boxplot(aes(fill=factor(Acc_Temp)))+
  facet_wrap(~day)+
  scale_fill_manual('Acclimated\nTemperature',values=c("blue","orange","red"))+
  theme_bw()

an.25.0<-glm(med~Acc_Temp,family=Gamma(link = "log"),weights = weights[bob.25$day==0],data=bob.25[bob.25$day==0,])
summary(an.25.0)

an.25.7<-glm(med~Acc_Temp,family=Gamma(link = "log"),weights = weights[bob.25$day==7],data=bob.25[bob.25$day==7,])
summary(an.25.7)

an.25.14<-glm(med~Acc_Temp,family=Gamma(link = "log"),weights = weights[bob.25$day==14],data=bob.25[bob.25$day==14,])
summary(an.25.14)



#### 37.5 C Acute: ####
bob.37<-res3[res3$Acute_Temp==37.5,]
#bob.37<-bob.37[!(bob.37$day==2 & bob.37$Acc_Temp==25),]
bob.37$Acc_Temp<-factor(bob.37$Acc_Temp,levels=c('37.5','12.5','25'))
bob.37<-na.omit(bob.37)

weights <- 1/bob.37$med.ci.width
weights <- weights/mean(weights)

gm.37.4<-gam(med~s(day,by=Acc_Temp,k=4)+Acc_Temp,family=Gamma(link="log"),data=bob.37,weights = weights)
gm.37.5<-gam(med~s(day,by=Acc_Temp,k=5)+Acc_Temp,family=Gamma(link="log"),data=bob.37,weights = weights)
gm.37.6<-gam(med~s(day,by=Acc_Temp,k=6)+Acc_Temp,family=Gamma(link="log"),data=bob.37,weights = weights)
gm.37.7<-gam(med~s(day,by=Acc_Temp,k=7)+Acc_Temp,family=Gamma(link="log"),data=bob.37,weights = weights)
gm.37.8<-gam(med~s(day,by=Acc_Temp,k=8)+Acc_Temp,family=Gamma(link="log"),data=bob.37,weights = weights)

# favors k = 7, but on visual inspection, badly over-fit
AICctab(gm.37.4,gm.37.5,gm.37.6,gm.37.7,gm.37.8)

gm.37<-gam(med~s(day,by=Acc_Temp,k=6)+Acc_Temp,family=Gamma(link="log"),data=bob.37,weights = weights)
summary(gm.37)

gm.37.0<-gam(med~s(day,k=6)+Acc_Temp,family=Gamma(link="log"),data=bob.37,weights = weights)
summary(gm.37.0)

# so, still support for distinct smooths
AICctab(gm.37,gm.37.0)



# Family: Gamma 
# Link function: log 
# 
# Formula:
#   med ~ s(day, by = Acc_Temp, k = 6) + Acc_Temp
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.391073   0.030386  45.779   <2e-16 ***
#   Acc_Temp12.5 0.097185   0.050843   1.911   0.0615 .  
# Acc_Temp25   0.007251   0.056670   0.128   0.8987    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df      F  p-value    
# s(day):Acc_Temp37.5 4.953  4.999  8.432 7.35e-06 ***
#   s(day):Acc_Temp12.5 4.784  4.975  3.974  0.00371 ** 
#   s(day):Acc_Temp25   4.630  4.929 15.174  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.734   Deviance explained = 76.1%
# GCV = 0.24687  Scale est. = 0.18653   n = 69


bob.37 %>% filter(day %in% c(0,7,14)) %>%
  ggplot(aes(y=med))+
  geom_boxplot(aes(fill=factor(Acc_Temp)))+
  facet_wrap(~day)+
  scale_fill_manual('Acclimated\nTemperature',values=c("blue","orange","red"))+
  theme_bw()

an.37.0<-glm(med~Acc_Temp,family=Gamma(link = "log"),weights = weights[bob.37$day==0],data=bob.37[bob.37$day==0,])
summary(an.37.0)

an.37.7<-glm(med~Acc_Temp,family=Gamma(link = "log"),weights = weights[bob.37$day==7],data=bob.37[bob.37$day==7,])
summary(an.37.7)

an.37.14<-glm(med~Acc_Temp,family=Gamma(link = "log"),weights = weights[bob.37$day==14],data=bob.37[bob.37$day==14,])
summary(an.37.14)





#### Visualizations: ####

threshold<-5*0.32029244

# Acute 12.5:

dat.12<-expand.grid(Acc_Temp=c('12.5','25','37.5'),Acute_Temp='12.5',
                    day=seq(0,14,0.2))
dat.12$Acc_Temp<-factor(dat.12$Acc_Temp,levels=c('12.5','25','37.5'))

pds12<-predict(gm.12,newdata = dat.12,type='response',se.fit = T)
dat.12$med<-pds12$fit
dat.12$med.se<-pds12$se.fit

p.12<-ggplot(dat.12,aes(x=day))+
  geom_line(aes(y=med,color=Acc_Temp))+
  geom_ribbon(aes(ymin=med-1.96*med.se,ymax=med+1.96*med.se,group=Acc_Temp),alpha=0.1)+
  geom_point(data=bob.12,aes(y=med,color=Acc_Temp),position=position_dodge(width=0.2))+
  geom_errorbar(data=bob.12,aes(ymin=med.lwr,ymax=med.upr,color=Acc_Temp),width=0.1,position = position_dodge(width=0.2))+
  geom_hline(yintercept = threshold,linetype=2)+
  scale_color_manual('Acclimated\nTemperature',values=c("blue","orange","red"))+
  scale_x_continuous('Day')+
  scale_y_continuous('Median velocity',limits=c(0,9),expand = c(0,0))+
  theme_bw()+
  theme(legend.position = c(0.87, 0.83),
        legend.background = element_rect(fill = "white", color = "black"))+
  ggtitle('Acute Temp = 12.5 C')
p.12


# Acute 25:

dat.25<-expand.grid(Acc_Temp=c('12.5','25','37.5'),Acute_Temp='25',
                    day=seq(0,14,0.2))
dat.25$Acc_Temp<-factor(dat.25$Acc_Temp,levels=c('25','12.5','37.5'))

pds25<-predict(gm.25,newdata = dat.25,type='response',se.fit = T)
dat.25$med<-pds25$fit
dat.25$med.se<-pds25$se.fit

p.25<-ggplot(dat.25,aes(x=day))+
  geom_line(aes(y=med,color=Acc_Temp))+
  geom_ribbon(aes(ymin=med-1.96*med.se,ymax=med+1.96*med.se,group=Acc_Temp),alpha=0.1)+
  geom_point(data=bob.25,aes(y=med,color=Acc_Temp),position=position_dodge(width=0.2))+
  geom_errorbar(data=bob.25,aes(ymin=med.lwr,ymax=med.upr,color=Acc_Temp),width=0.1,position = position_dodge(width=0.2))+
  geom_hline(yintercept = threshold,linetype=2)+
  scale_color_manual('Acclimated\nTemperature',values=c("orange","blue","red"))+
  scale_x_continuous('Day')+
  scale_y_continuous('Median velocity',limits=c(0,9),expand = c(0,0))+
  theme_bw()+
#  theme(legend.position = c(0.87, 0.85),
#        legend.background = element_rect(fill = "white", color = "black"))+
  theme(legend.position = "none")+
  ggtitle('Acute Temp = 25 C')
p.25


# Acute 37.5:

dat.37<-expand.grid(Acc_Temp=c('12.5','25','37.5'),Acute_Temp='37.5',
                    day=seq(0,14,0.2))
dat.37$Acc_Temp<-factor(dat.37$Acc_Temp,levels=c('37.5','12.5','25'))

pds37<-predict(gm.37,newdata = dat.37,type='response',se.fit = T)
dat.37$med<-pds37$fit
dat.37$med.se<-pds37$se.fit

p.37<-ggplot(dat.37,aes(x=day))+
  geom_line(aes(y=med,color=Acc_Temp))+
  geom_ribbon(aes(ymin=med-1.96*med.se,ymax=med+1.96*med.se,group=Acc_Temp),alpha=0.1)+
  geom_point(data=bob.37,aes(y=med,color=Acc_Temp),position=position_dodge(width=0.2))+
  geom_errorbar(data=bob.37,aes(ymin=med.lwr,ymax=med.upr,color=Acc_Temp),width=0.1,position = position_dodge(width=0.2))+
  geom_hline(yintercept = threshold,linetype=2)+
  scale_color_manual('Acclimated\nTemperature',values=c("red","blue","orange"))+
  scale_x_continuous('Day')+
  scale_y_continuous('Median velocity',limits=c(0,9),expand = c(0,0))+
  theme_bw()+
#  theme(legend.position = c(0.87, 0.85),
#        legend.background = element_rect(fill = "white", color = "black"))+
  theme(legend.position = "none")+
  ggtitle('Acute Temp = 37.5 C')
p.37

grid.arrange(p.12,p.25,p.37,nrow=1)

#### Save final plot: ####

a1<-arrangeGrob(p.12,p.25,p.37,nrow=1)
scl<-1.2
ggsave(plot = a1,filename = "./results/Median_speed_gams_v2.pdf",width=scl*13,height=scl*4)




#### B. Proportion moving analysis #####

# Round data up to avoid true 0's that are incompatible with beta distributions:
#min(res3$pmove[res3$pmove!=0])
# 0.0006612028
res3$pmove<-ifelse(res3$pmove==0,0.0001,res3$pmove)


#### Difference in overall activity by Acute temp for fully acclimated populations ###

tmp<-res3 %>% filter(Acute_Temp==Acc_Temp)

j0<-gam(pmove~1,family=betar,data=tmp)
j1<-gam(pmove~factor(Acc_Temp),family=betar,data=tmp)
summary(j1)

AICctab(j0,j1)

expit(coef(j1)[1])
expit(coef(j1)[1]+coef(j1)[2])
expit(coef(j1)[1]+coef(j1)[3])



#### 12 C Acute: ####
bob.12<-res3[res3$Acute_Temp==12.5,]
bob.12$Acc_Temp<-as.factor(bob.12$Acc_Temp)

weights <- 1/bob.12$pmove.ci.width
weights <- weights/mean(weights)

gm.12.4<-gam(pmove~s(day,by=Acc_Temp,k=4)+Acc_Temp,family=betar,data=bob.12,weights = weights)
gm.12.5<-gam(pmove~s(day,by=Acc_Temp,k=5)+Acc_Temp,family=betar,data=bob.12,weights = weights)
gm.12.6<-gam(pmove~s(day,by=Acc_Temp,k=6)+Acc_Temp,family=betar,data=bob.12,weights = weights)
gm.12.7<-gam(pmove~s(day,by=Acc_Temp,k=7)+Acc_Temp,family=betar,data=bob.12,weights = weights)
gm.12.8<-gam(pmove~s(day,by=Acc_Temp,k=8)+Acc_Temp,family=betar,data=bob.12,weights = weights)

# favors all the k's have similar support, will go with 6 for consistency
AICctab(gm.12.4,gm.12.5,gm.12.6,gm.12.7,gm.12.8)

# use k=6
gm.12<-gam(pmove~s(day,by=Acc_Temp,k=6)+Acc_Temp,family=betar,data=bob.12,weights = weights)
summary(gm.12)

gm.12.0<-gam(pmove~s(day,k=6)+Acc_Temp,family=betar,data=bob.12,weights = weights)
summary(gm.12.0)

# so, distinct smooths are not strongly supported here!
AICctab(gm.12,gm.12.0)


#### 25 C Acute: ####
bob.25<-res3[res3$Acute_Temp==25,]
bob.25$Acc_Temp<-factor(bob.25$Acc_Temp,levels=c('25','12.5','37.5'))

weights <- 1/bob.25$pmove.ci.width
weights <- weights/mean(weights)

gm.25.4<-gam(pmove~s(day,by=Acc_Temp,k=4)+Acc_Temp,family=betar,data=bob.25,weights = weights)
gm.25.5<-gam(pmove~s(day,by=Acc_Temp,k=5)+Acc_Temp,family=betar,data=bob.25,weights = weights)
gm.25.6<-gam(pmove~s(day,by=Acc_Temp,k=6)+Acc_Temp,family=betar,data=bob.25,weights = weights)
gm.25.7<-gam(pmove~s(day,by=Acc_Temp,k=7)+Acc_Temp,family=betar,data=bob.25,weights = weights)
gm.25.8<-gam(pmove~s(day,by=Acc_Temp,k=8)+Acc_Temp,family=betar,data=bob.25,weights = weights)

# favors k = 7, but on visual inspection, over-fit
AICctab(gm.25.4,gm.25.5,gm.25.6,gm.25.7,gm.25.8)

# use k=6
gm.25<-gam(pmove~s(day,by=Acc_Temp,k=6)+Acc_Temp,family=betar,data=bob.25,weights = weights)
summary(gm.25)

gm.25.0<-gam(pmove~s(day,k=6)+Acc_Temp,family=betar,data=bob.25,weights = weights)
summary(gm.25.0)

# this one has support for distinct smooths
AICctab(gm.25,gm.25.0)


#### 37.5 C Acute: ####
bob.37<-res3[res3$Acute_Temp==37.5,]
#bob.37<-bob.37[!(bob.37$day==2 & bob.37$Acc_Temp==25),]
bob.37$Acc_Temp<-factor(bob.37$Acc_Temp,levels=c('37.5','12.5','25'))
bob.37<-na.omit(bob.37)
nrow(bob.37)

weights <- 1/bob.37$pmove.ci.width
weights <- weights/mean(weights)

gm.37.4<-gam(pmove~s(day,by=Acc_Temp,k=4)+Acc_Temp,family=betar,data=bob.37,weights = weights)
gm.37.5<-gam(pmove~s(day,by=Acc_Temp,k=5)+Acc_Temp,family=betar,data=bob.37,weights = weights)
gm.37.6<-gam(pmove~s(day,by=Acc_Temp,k=6)+Acc_Temp,family=betar,data=bob.37,weights = weights)
gm.37.7<-gam(pmove~s(day,by=Acc_Temp,k=7)+Acc_Temp,family=betar,data=bob.37,weights = weights)
gm.37.8<-gam(pmove~s(day,by=Acc_Temp,k=8)+Acc_Temp,family=betar,data=bob.37,weights = weights)

# favors k = 5, but not strongly; sticking with k=6
AICctab(gm.37.4,gm.37.5,gm.37.6,gm.37.7,gm.37.8)

gm.37<-gam(pmove~s(day,by=Acc_Temp,k=6)+Acc_Temp,family=betar,data=bob.37,weights = weights)
summary(gm.37)

gm.37.0<-gam(pmove~s(day,k=6)+Acc_Temp,family=betar,data=bob.37,weights = weights)
summary(gm.37.0)

# so, moderate support for distinct smooths
AICctab(gm.37,gm.37.0)


#### Visualizations: ####

# Acute 12.5:
dat.12<-expand.grid(Acc_Temp=c('12.5','25','37.5'),Acute_Temp='12.5',
                    day=seq(0,14,0.2))
dat.12$Acc_Temp<-factor(dat.12$Acc_Temp,levels=c('12.5','25','37.5'))

pds12<-predict(gm.12,newdata = dat.12,type='link',se.fit = T)
dat.12$pmove<-pds12$fit
dat.12$pmove.lwr<-pds12$fit-1.96*pds12$se.fit
dat.12$pmove.upr<-pds12$fit+1.96*pds12$se.fit

dat.12$pmove<-expit(dat.12$pmove)
dat.12$pmove.lwr<-expit(dat.12$pmove.lwr)
dat.12$pmove.upr<-expit(dat.12$pmove.upr)


p.12<-ggplot(dat.12,aes(x=day))+
  geom_line(aes(y=pmove,color=Acc_Temp))+
  geom_ribbon(aes(ymin=pmove.lwr,ymax=pmove.upr,group=Acc_Temp),alpha=0.1)+
  geom_point(data=bob.12,aes(y=pmove,color=Acc_Temp),position=position_dodge(width=0.2))+
  geom_errorbar(data=bob.12,aes(ymin=pmove.lwr,ymax=pmove.upr,color=Acc_Temp),width=0.1,position = position_dodge(width=0.2))+
  scale_color_manual('Acclimated\nTemperature',values=c("blue","orange","red"))+
  scale_x_continuous('Day',limits=c(-0.1,15))+
  scale_y_continuous('Proportion moving',limits=c(0,1))+
  theme_bw()+
  theme(legend.position = c(0.87, 0.83),
        legend.background = element_rect(fill = "white", color = "black"))+
  ggtitle('Acute Temp = 12.5 C')
p.12


# Acute 25:
dat.25<-expand.grid(Acc_Temp=c('12.5','25','37.5'),Acute_Temp='25',
                    day=seq(0,14,0.2))
dat.25$Acc_Temp<-factor(dat.25$Acc_Temp,levels=c('25','12.5','37.5'))

pds25<-predict(gm.25,newdata = dat.25,type='link',se.fit = T)
dat.25$pmove<-pds25$fit
dat.25$pmove.lwr<-pds25$fit-1.96*pds25$se.fit
dat.25$pmove.upr<-pds25$fit+1.96*pds25$se.fit

dat.25$pmove<-expit(dat.25$pmove)
dat.25$pmove.lwr<-expit(dat.25$pmove.lwr)
dat.25$pmove.upr<-expit(dat.25$pmove.upr)


p.25<-ggplot(dat.25,aes(x=day))+
  geom_line(aes(y=pmove,color=Acc_Temp))+
  geom_ribbon(aes(ymin=pmove.lwr,ymax=pmove.upr,group=Acc_Temp),alpha=0.1)+
  geom_point(data=bob.25,aes(y=pmove,color=Acc_Temp),position=position_dodge(width=0.2))+
  geom_errorbar(data=bob.25,aes(ymin=pmove.lwr,ymax=pmove.upr,color=Acc_Temp),width=0.1,position = position_dodge(width=0.2))+
  scale_color_manual('Acclimated\nTemperature',values=c("orange","blue","red"))+
  scale_x_continuous('Day',limits=c(-0.1,15))+
  scale_y_continuous('Proportion moving',limits=c(0,1))+
  theme_bw()+
  #  theme(legend.position = c(0.87, 0.85),
  #        legend.background = element_rect(fill = "white", color = "black"))+
  theme(legend.position = "none")+
  ggtitle('Acute Temp = 25 C')
p.25


# Acute 37.5:
dat.37<-expand.grid(Acc_Temp=c('12.5','25','37.5'),Acute_Temp='37.5',
                    day=seq(0,14,0.2))
dat.37$Acc_Temp<-factor(dat.37$Acc_Temp,levels=c('37.5','12.5','25'))

pds37<-predict(gm.37,newdata = dat.37,type='link',se.fit = T)
dat.37$pmove<-pds37$fit
dat.37$pmove.lwr<-pds37$fit-1.96*pds37$se.fit
dat.37$pmove.upr<-pds37$fit+1.96*pds37$se.fit

dat.37$pmove<-expit(dat.37$pmove)
dat.37$pmove.lwr<-expit(dat.37$pmove.lwr)
dat.37$pmove.upr<-expit(dat.37$pmove.upr)


p.37<-ggplot(dat.37,aes(x=day))+
  geom_line(aes(y=pmove,color=Acc_Temp))+
  geom_ribbon(aes(ymin=pmove.lwr,ymax=pmove.upr,group=Acc_Temp),alpha=0.1)+
  geom_point(data=bob.37,aes(y=pmove,color=Acc_Temp),position=position_dodge(width=0.2))+
  geom_errorbar(data=bob.37,aes(ymin=pmove.lwr,ymax=pmove.upr,color=Acc_Temp),width=0.1,position = position_dodge(width=0.2))+
  scale_color_manual('Acclimated\nTemperature',values=c("red","blue","orange"))+
  scale_x_continuous('Day',limits=c(-0.1,15))+
  scale_y_continuous('Proportion moving',limits=c(0,1))+
  theme_bw()+
  #  theme(legend.position = c(0.87, 0.85),
  #        legend.background = element_rect(fill = "white", color = "black"))+
  theme(legend.position = "none")+
  ggtitle('Acute Temp = 37.5 C')
p.37


grid.arrange(p.12,p.25,p.37,nrow=1)

#### Save final plot: ####
a1<-arrangeGrob(p.12,p.25,p.37,nrow=1)
scl<-1.2
ggsave(plot = a1,filename = "./results/pmove_gams_v2.pdf",width=scl*13,height=scl*4)




#### C. Turning angle ####


#### 12 C Acute: ####
bob.12<-res3[res3$Acute_Temp==12.5,]
bob.12$Acc_Temp<-as.factor(bob.12$Acc_Temp)

weights <- 1/bob.12$mu.ci.width
weights <- weights/mean(weights)

gm.12.4<-gam(mu~s(day,by=Acc_Temp,k=4)+Acc_Temp,family=betar,data=bob.12,weights = weights)
gm.12.5<-gam(mu~s(day,by=Acc_Temp,k=5)+Acc_Temp,family=betar,data=bob.12,weights = weights)
gm.12.6<-gam(mu~s(day,by=Acc_Temp,k=6)+Acc_Temp,family=betar,data=bob.12,weights = weights)
gm.12.7<-gam(mu~s(day,by=Acc_Temp,k=7)+Acc_Temp,family=betar,data=bob.12,weights = weights)
gm.12.8<-gam(mu~s(day,by=Acc_Temp,k=8)+Acc_Temp,family=betar,data=bob.12,weights = weights)

# most of the k's have similar support, will go with 6 for consistency
AICctab(gm.12.4,gm.12.5,gm.12.6,gm.12.7,gm.12.8)

# use k=6
gm.12<-gam(mu~s(day,by=Acc_Temp,k=6)+Acc_Temp,family=betar,data=bob.12,weights = weights)
summary(gm.12)

gm.12.0<-gam(mu~s(day,k=6)+Acc_Temp,family=betar,data=bob.12,weights = weights)
summary(gm.12.0)

# so, distinct smooths are not strongly supported here!
AICctab(gm.12,gm.12.0)


#### 25 C Acute: ####
bob.25<-res3[res3$Acute_Temp==25,]
bob.25$Acc_Temp<-factor(bob.25$Acc_Temp,levels=c('25','12.5','37.5'))

weights <- 1/bob.25$mu.ci.width
weights <- weights/mean(weights)

gm.25.4<-gam(mu~s(day,by=Acc_Temp,k=4)+Acc_Temp,family=betar,data=bob.25,weights = weights)
gm.25.5<-gam(mu~s(day,by=Acc_Temp,k=5)+Acc_Temp,family=betar,data=bob.25,weights = weights)
gm.25.6<-gam(mu~s(day,by=Acc_Temp,k=6)+Acc_Temp,family=betar,data=bob.25,weights = weights)
gm.25.7<-gam(mu~s(day,by=Acc_Temp,k=7)+Acc_Temp,family=betar,data=bob.25,weights = weights)
gm.25.8<-gam(mu~s(day,by=Acc_Temp,k=8)+Acc_Temp,family=betar,data=bob.25,weights = weights)

# favors k = 8, but on visual inspection, over-fit
AICctab(gm.25.4,gm.25.5,gm.25.6,gm.25.7,gm.25.8)

# use k=6
gm.25<-gam(mu~s(day,by=Acc_Temp,k=6)+Acc_Temp,family=betar,data=bob.25,weights = weights)
summary(gm.25)

gm.25.0<-gam(mu~s(day,k=6)+Acc_Temp,family=betar,data=bob.25,weights = weights)
summary(gm.25.0)

# this one has support for distinct smooths
AICctab(gm.25,gm.25.0)


#### 37.5 C Acute: ####
bob.37<-res3[res3$Acute_Temp==37.5,]
#bob.37<-bob.37[!(bob.37$day==2 & bob.37$Acc_Temp==25),]
bob.37$Acc_Temp<-factor(bob.37$Acc_Temp,levels=c('37.5','12.5','25'))
bob.37<-na.omit(bob.37)
nrow(bob.37)

weights <- 1/bob.37$mu.ci.width
weights <- weights/mean(weights)

gm.37.4<-gam(mu~s(day,by=Acc_Temp,k=4)+Acc_Temp,family=betar,data=bob.37,weights = weights)
gm.37.5<-gam(mu~s(day,by=Acc_Temp,k=5)+Acc_Temp,family=betar,data=bob.37,weights = weights)
gm.37.6<-gam(mu~s(day,by=Acc_Temp,k=6)+Acc_Temp,family=betar,data=bob.37,weights = weights)
gm.37.7<-gam(mu~s(day,by=Acc_Temp,k=7)+Acc_Temp,family=betar,data=bob.37,weights = weights)
gm.37.8<-gam(mu~s(day,by=Acc_Temp,k=8)+Acc_Temp,family=betar,data=bob.37,weights = weights)

# favors k = 5, but not strongly; sticking with k=6
AICctab(gm.37.4,gm.37.5,gm.37.6,gm.37.7,gm.37.8)

gm.37<-gam(mu~s(day,by=Acc_Temp,k=6)+Acc_Temp,family=betar,data=bob.37,weights = weights)
summary(gm.37)

gm.37.0<-gam(mu~s(day,k=6)+Acc_Temp,family=betar,data=bob.37,weights = weights)
summary(gm.37.0)

# so, moderate support for distinct smooths
AICctab(gm.37,gm.37.0)


#### Visualizations: ####

# Acute 12.5:
dat.12<-expand.grid(Acc_Temp=c('12.5','25','37.5'),Acute_Temp='12.5',
                    day=seq(0,14,0.2))
dat.12$Acc_Temp<-factor(dat.12$Acc_Temp,levels=c('12.5','25','37.5'))

pds12<-predict(gm.12,newdata = dat.12,type='link',se.fit = T)
dat.12$mu<-pds12$fit
dat.12$mu.lwr<-pds12$fit-1.96*pds12$se.fit
dat.12$mu.upr<-pds12$fit+1.96*pds12$se.fit

dat.12$mu<-expit(dat.12$mu)
dat.12$mu.lwr<-expit(dat.12$mu.lwr)
dat.12$mu.upr<-expit(dat.12$mu.upr)


p.12<-ggplot(dat.12,aes(x=day))+
  geom_line(aes(y=mu,color=Acc_Temp))+
  geom_ribbon(aes(ymin=mu.lwr,ymax=mu.upr,group=Acc_Temp),alpha=0.1)+
  geom_point(data=bob.12,aes(y=mu,color=Acc_Temp),position=position_dodge(width=0.2))+
  geom_errorbar(data=bob.12,aes(ymin=mu.lwr,ymax=mu.upr,color=Acc_Temp),width=0.1,position = position_dodge(width=0.2))+
  scale_color_manual('Acclimated\nTemperature',values=c("blue","orange","red"))+
  scale_x_continuous('Day',limits=c(-0.1,15))+
  scale_y_continuous('Scaled absolute angle',limits=c(0,1))+
  theme_bw()+
  theme(legend.position = c(0.87, 0.83),
        legend.background = element_rect(fill = "white", color = "black"))+
  ggtitle('Acute Temp = 12.5 C')
p.12


# Acute 25:
dat.25<-expand.grid(Acc_Temp=c('12.5','25','37.5'),Acute_Temp='25',
                    day=seq(0,14,0.2))
dat.25$Acc_Temp<-factor(dat.25$Acc_Temp,levels=c('25','12.5','37.5'))

pds25<-predict(gm.25,newdata = dat.25,type='link',se.fit = T)
dat.25$mu<-pds25$fit
dat.25$mu.lwr<-pds25$fit-1.96*pds25$se.fit
dat.25$mu.upr<-pds25$fit+1.96*pds25$se.fit

dat.25$mu<-expit(dat.25$mu)
dat.25$mu.lwr<-expit(dat.25$mu.lwr)
dat.25$mu.upr<-expit(dat.25$mu.upr)


p.25<-ggplot(dat.25,aes(x=day))+
  geom_line(aes(y=mu,color=Acc_Temp))+
  geom_ribbon(aes(ymin=mu.lwr,ymax=mu.upr,group=Acc_Temp),alpha=0.1)+
  geom_point(data=bob.25,aes(y=mu,color=Acc_Temp),position=position_dodge(width=0.2))+
  geom_errorbar(data=bob.25,aes(ymin=mu.lwr,ymax=mu.upr,color=Acc_Temp),width=0.1,position = position_dodge(width=0.2))+
  scale_color_manual('Acclimated\nTemperature',values=c("orange","blue","red"))+
  scale_x_continuous('Day',limits=c(-0.1,15))+
  scale_y_continuous('Scaled absolute angle',limits=c(0,1))+
  theme_bw()+
  #  theme(legend.position = c(0.87, 0.85),
  #        legend.background = element_rect(fill = "white", color = "black"))+
  theme(legend.position = "none")+
  ggtitle('Acute Temp = 25 C')
p.25


# Acute 37.5:
dat.37<-expand.grid(Acc_Temp=c('12.5','25','37.5'),Acute_Temp='37.5',
                    day=seq(0,14,0.2))
dat.37$Acc_Temp<-factor(dat.37$Acc_Temp,levels=c('37.5','12.5','25'))

pds37<-predict(gm.37,newdata = dat.37,type='link',se.fit = T)
dat.37$mu<-pds37$fit
dat.37$mu.lwr<-pds37$fit-1.96*pds37$se.fit
dat.37$mu.upr<-pds37$fit+1.96*pds37$se.fit

dat.37$mu<-expit(dat.37$mu)
dat.37$mu.lwr<-expit(dat.37$mu.lwr)
dat.37$mu.upr<-expit(dat.37$mu.upr)


p.37<-ggplot(dat.37,aes(x=day))+
  geom_line(aes(y=mu,color=Acc_Temp))+
  geom_ribbon(aes(ymin=mu.lwr,ymax=mu.upr,group=Acc_Temp),alpha=0.1)+
  geom_point(data=bob.37,aes(y=mu,color=Acc_Temp),position=position_dodge(width=0.2))+
  geom_errorbar(data=bob.37,aes(ymin=mu.lwr,ymax=mu.upr,color=Acc_Temp),width=0.1,position = position_dodge(width=0.2))+
  scale_color_manual('Acclimated\nTemperature',values=c("red","blue","orange"))+
  scale_x_continuous('Day',limits=c(-0.1,15))+
  scale_y_continuous('Scaled absolute angle',limits=c(0,1))+
  theme_bw()+
  #  theme(legend.position = c(0.87, 0.85),
  #        legend.background = element_rect(fill = "white", color = "black"))+
  theme(legend.position = "none")+
  ggtitle('Acute Temp = 37.5 C')
p.37


grid.arrange(p.12,p.25,p.37,nrow=1)

#### Save final plot: ####
a1<-arrangeGrob(p.12,p.25,p.37,nrow=1)
scl<-1.2
ggsave(plot = a1,filename = "./results/mu_gams_v2.pdf",width=scl*13,height=scl*4)


