###Simulation study
#NOTE 1: in Rstudio; press ALT+o for overview of code section headings

#Load libraries----
setwd("C:\\Users\\VNOB-0728\\Documents\\Server\\Analysis squeezy")
pathResults <- "./SimulationResults/" 
pathFigures <- "./SimulationFigures/" 
#setwd("C:\\Synchr\\Rscripts\\EB\\Squeezy")
library(ecpc)
library(MASS)
library(penalized)
library(glmnet)
library(mvtnorm)
library(mgcv)
library(CVXR)
library(foreach)
library(doParallel)

library(gren)
library(ipflasso)
library(fwelnet)

library(multiridge)
source("multiridge.R")
source("squeezy.R")

#library(Matrix)
#library(GRridge)

#load libraries needed for storing and plotting results
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(scales)
library(tidyr)

#optional: setup parallel backend to use many processors
# cores=detectCores()
if(0){
  cl <- makeCluster(6) #not to overload your computer
  registerDoParallel(cl)
}

#generate seeds for reproducibility
#set.seed(1)
seeds <- round(runif(10000)*10^6)
#save(seeds,file=paste(pathResults,"seeds.R",sep=""))
load(paste(pathResults,"seeds.R",sep=""))

#Set simulation variables----
nSim <- 100
n<-150
p<-1200
n2<-300 #sample size of each test data set
q<-0 #set q*100% of the betas to 0
alphaEN <- 0.3

#Plots: general parameters----
wdth<-600
hght<-wdth*5/8
wdthpdf <- wdth/75
hghtpdf <- hght/75
ts <- 16 #basis text size in figures
ls <- 1.5 #basis line size in figures
ps <- 2 #basis point size in figures
sz <- 2 #point size
strk <- 1.5 #stroke size
palette <- "Dark2"
#display.brewer.all(3,colorblindFriendly=T)
colpal <- "Dark2"
colsfill <- brewer.pal(3,"Dark2")[1:2]
colsAUC <- brewer.pal(8,"Dark2")
colsMSE <- seq_gradient_pal(low="grey50",high="white")((0:8)/9)
colsTime <- seq_gradient_pal(low="black",high="grey50")((0:8)/9)


#load data----
# fname <- paste(pathResults,paste("simres",alphaEN,p,n,taugrp[1],taugrp[5],sigmasq,sep="_"),".Rdata",sep="")
# load(fname)

#plots performance MSE
alphaEN <- 0.3
df2 <- data.frame()
for(alphaEN in c(0.3,0.8)){
  for(setting in c("1","2","3")){
    fname <- paste(pathResults,paste("simres",alphaEN,setting,sep="_"),".Rdata",sep="")
    #if(setting=="1"&alphaEN==0.3) fname <- paste(pathResults,paste("simreshat",alphaEN,setting,sep="_"),".Rdata",sep="")
    load(fname)
    #df <- df[!grepl("sq",df$Method),]
    #df <- df[!df$Method=="ecpc",]
    df$Setting <- setting
    df$Alpha <- alphaEN
    df2 <- rbind(df2,df)
    
    fname <- paste(pathResults,paste("simresSqueezy",alphaEN,setting,sep="_"),".Rdata",sep="")
    load(fname)
    df$Setting <- setting
    df$Alpha <- alphaEN
    df2 <- rbind(df2,df)
  }
}

df2$Method <- factor(df2$Method,levels=unique(df2$Method)[c(1,3,2,8,6,4,7,5)],
                        labels=c("EN","fwEN","ipf","ecpcEN squeezy","squeezy (single)",
                                 "squeezy (multi)","squeezy (single+reCV)","squeezy (multi+reCV)"))

df2$Setting2 <- as.numeric(df2$Setting)
df2$Setting <- factor(df2$Setting,levels=c("1","2","3"),
                      labels=c("No groups","Weakly informative groups","Informative groups"))
df2$MethodNum <- as.factor(as.numeric(df2$Method))
df2$MethodBoth <- as.factor(paste(df2$MethodNum,". ",df2$Method,sep=""))

summBox <- df2 %>% group_by(Setting2,Alpha) %>% summarise("low"=boxplot.stats(MSE)$stats[1],
                                                    "high"=boxplot.stats(MSE)$stats[5]) %>% ungroup()
df2 <- df2 %>% left_join(summBox,by=c("Alpha","Setting2"))
df2NoOutliers <- df2[df2$MSE>=df2$low & df2$MSE<=df2$high,]


#plots time
df3 <- data.frame()
alphaEN <- 0.3
for(setting in as.character(4:11)){
  fname <- paste(pathResults,paste("simres",alphaEN,setting,sep="_"),".Rdata",sep="")
  load(fname)
  df <- df[!grepl("sq",df$Method),]
  df <- df[!df$Method=="ecpc",]
  df$Setting <- setting
  if(setting=="3"){
    df <- df[df$Dataset%in%1:5,]
  }
  df3 <- rbind(df3,df)
  
  fname <- paste(pathResults,paste("simresSqueezy",alphaEN,setting,sep="_"),".Rdata",sep="")
  load(fname)
  df$Setting <- setting
  df3 <- rbind(df3,df)
  
  if(setting%in%c(4,6,8,10)){
    fname <- paste(pathResults,paste("simresIPF2",alphaEN,setting,sep="_"),".Rdata",sep="")
    load(fname)
    df$Setting <- setting
    df3 <- rbind(df3,df)
  }
}
df3 <- df3[!is.nan(df3$time),]
df3$Method <- factor(df3$Method,levels=unique(df3$Method)[c(1,3,2,9,8,6,4,7,5)],
                     labels=c("EN","fwEN","ipf","ipf2","ecpcEN squeezy",
                              "squeezy (single)","squeezy (multi)",
                              "squeezy (single+reCV)","squeezy (multi+reCV)"))
df3$MethodNum <- as.factor(as.numeric(as.factor(df3$Method)))
df3$MethodBoth <- as.factor(paste(df3$MethodNum,". ",df3$Method,sep=""))

summdf3 <- df3 %>% group_by(Setting,Method,n,G,p,MethodBoth,MethodNum) %>% summarise(meanTime=mean(time)) %>% ungroup()


#plot alpha lambda
alphaEN <- 0.3
setting <- "3"
fname <- paste(pathResults,paste("simressqueezyalpha",setting,sep="_"),".Rdata",sep="")
load(fname) #df, dfgrps

df4 <- df

Summdf4 <- df4 %>% group_by(reCV,Alpha) %>% summarise(mean=mean(MSEApprox),
                                                      median=median(MSEApprox),
                                                      q25=quantile(MSEApprox,0.25),q75=quantile(MSEApprox,0.75),
                                                      q95=quantile(MSEApprox,0.95),q05=quantile(MSEApprox,0.05)) %>% ungroup()
Summdf4$reCV <- factor(Summdf4$reCV,levels = c(TRUE,FALSE),labels=c("TRUE","FALSE"))

#Plot: prediction all with outliers----
alphaEN <- 0.3
figname<-paste(pathFigures,"FigSimsMSE",alphaEN,"AllOutliers.pdf",sep="")
pdf(width = wdthpdf*1.5, height = hghtpdf*1.1,
    file = figname)
p1<-ggplot(df2[df2$Alpha==alphaEN,])+aes(x=MethodNum,y=MSE)+
  geom_boxplot(aes(fill=MethodBoth))+
  facet_wrap(.~Setting,scales="free")+
  scale_fill_manual(values=colsMSE,name="Method")+
  labs(y="MSE",x="Method")+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts),
        legend.position ="bottom")#,
p1
dev.off()
p1
  

#Plot: prediction performance all no outliers----
alphaEN <- 0.3
figname<-paste(pathFigures,"FigSimsMSE",alphaEN,"AllNoOutliers.pdf",sep="")
pdf(width = wdthpdf*1.7, height = hghtpdf,
    file = figname)
p1<-ggplot(df2NoOutliers[df2NoOutliers$Alpha==alphaEN,])+aes(x=MethodNum,y=MSE)+
  geom_boxplot(aes(fill=MethodBoth),outlier.shape = NA)+
  facet_wrap(.~Setting,scales="free")+
  labs(y="MSE",x="Method")+
  scale_fill_manual(values=colsMSE,name="Method")+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts),
        legend.position ="bottom")#,
p1
dev.off()
p1


#Plot: time dependence n----
figname<-paste(pathFigures,"FigSimsTimeN.pdf",sep="")
pdf(width = wdthpdf, height = hghtpdf,
    file = figname)
p1<-ggplot(df3[df3$Setting%in%c(10,4,5),])+
  aes(x=n,y=time,col=MethodBoth)+
  #geom_jitter(width=5,height=0,alpha=0.8,size=ps)+
  geom_point(data=summdf3[summdf3$Setting%in%c(10,4,5),],
             aes(y=meanTime,col=MethodBoth),size=ps)+
  geom_text_repel(data=summdf3[summdf3$Setting%in%c(10,4,5),],
            aes(y=meanTime,label=MethodNum),min.segment.length = 0.1,size=ts/3)+
  geom_line(data=summdf3[summdf3$Setting%in%c(10,4,5),],
            aes(y=meanTime),linetype=2,size=ls,alpha=0.2)+
  scale_color_manual(values=colsTime[1:9],name="Method")+
  #facet_wrap(.~n,scales="free_x")+
  labs(y="Time (s)",x="n")+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
p1
dev.off()
p1

#Plot: time dependence p----
figname<-paste(pathFigures,"FigSimsTimeP.pdf",sep="")
pdf(width = wdthpdf, height = hghtpdf,
    file = figname)
p2<-ggplot(df3[df3$Setting%in%c(10,6,7,11),])+
  aes(x=p,y=time,col=MethodBoth)+
  #geom_jitter(width=5,height=0,alpha=0.8,size=ps)+
  #geom_point()+
  geom_point(data=summdf3[summdf3$Setting%in%c(10,6,7,11),],
             aes(y=meanTime,col=MethodBoth),size=ps)+
  geom_text_repel(data=summdf3[summdf3$Setting%in%c(10,6,7,11),],
                  aes(y=meanTime,label=MethodNum),min.segment.length = 0.1,size=ts/3)+
  geom_line(data=summdf3[summdf3$Setting%in%c(10,6,7,11),],
            aes(y=meanTime),linetype=2,size=ls,alpha=0.2)+
  scale_color_manual(values=colsTime[1:9],name="Method")+
  scale_x_continuous(breaks = (1:4)*1000)+
  #facet_wrap(.~n,scales="free_x")+
  labs(y="")+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
p2
dev.off()
p2

#Plot: time dependence G----
figname<-paste(pathFigures,"FigSimsTimeG.pdf",sep="")
pdf(width = wdthpdf, height = hghtpdf,
    file = figname)
p3<-ggplot(df3[df3$Setting%in%c(10,8,9),])+
  aes(x=G,y=time,col=MethodBoth)+
  #geom_jitter(width=5,height=0,alpha=0.8,size=ps)+
  #geom_point()+
  geom_point(data=summdf3[summdf3$Setting%in%c(10,8,9),],
             aes(y=meanTime,col=MethodBoth),size=ps)+
  geom_text_repel(data=summdf3[summdf3$Setting%in%c(10,8,9),],
                  aes(y=meanTime,label=MethodNum),min.segment.length = 0.1,size=ts/3)+
  geom_line(data=summdf3[summdf3$Setting%in%c(10,8,9),],
            aes(y=meanTime),linetype=2,size=ls,alpha=0.2)+
  scale_color_manual(values=colsTime[1:9],name="Method")+
  #facet_wrap(.~n,scales="free_x")+
  labs(y="")+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
p3
dev.off()
p3

#Plot: time dependence n, p and G----
figname<-paste(pathFigures,"FigSimsTimeAll.pdf",sep="")
pdf(width = wdthpdf*1.7, height = hghtpdf,
    file = figname)
p4<-ggarrange(p1,p2,p3,common.legend=T,nrow=1,legend="bottom")+
  #theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
p4
dev.off()
p4
ggsave(figname,width= wdthpdf*1.7, height = hghtpdf)

#Plot: MSE different alpha----
figname<-paste(pathFigures,"FigMSEalpha.pdf",sep="")
pdf(width = wdthpdf, height = hghtpdf,
    file = figname)
p1<-ggplot(Summdf4)+
  aes(x=Alpha)+
  #geom_jitter(width=5,height=0,alpha=0.8,size=ps)+
  #geom_boxplot(aes(group=Alpha),fill="grey70")+
  geom_ribbon(aes(x=Alpha,ymin=q25,ymax=q75,fill=reCV),alpha=0.2)+
  #geom_ribbon(aes(x=Alpha,ymin=q05,ymax=q95,fill=reCV),alpha=0.1)+
  geom_line(aes(x=Alpha,y=q75,col=reCV,linetype=reCV),size=ls)+
  geom_line(aes(x=Alpha,y=q25,col=reCV,linetype=reCV),size=ls)+
  geom_line(aes(x=Alpha,y=q05,col=reCV,linetype=reCV),size=ls)+
  geom_line(aes(x=Alpha,y=q95,col=reCV,linetype=reCV),size=ls)+
  geom_line(aes(x=Alpha,y=median,col=reCV,linetype=reCV),size=ls)+
  scale_color_manual(values=c(colsMSE[2],"black"))+
  scale_fill_manual(values=c(colsMSE[2],"black"))+
  scale_linetype_manual(values=c(3,1))+
  # geom_line(data=df4[df4$Dataset%in%c(2),],
  #           aes(group=Dataset),
  #           size=ls,alpha=0.7,linetype=2)+
  #geom_smooth(col="black",linetype=3)+
  labs(y="MSE")+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
p1
dev.off()
p1

#Plot: lambda different alpha----
dfGrps$loglam <- log(dfGrps$lambdaApprox)
dfGrps$Group <- as.factor(dfGrps$Group)
dfGrps$GroupAlpha <- as.factor(paste(dfGrps$Alpha,dfGrps$Group,sep="."))
i<-2

SummLam <- dfGrps %>% group_by(Alpha,Group) %>% summarise(medianLogLam = median(loglam),
                                                          meanLogLam = mean(loglam),
                                                          q95=quantile(loglam,0.95),
                                                          q05=quantile(loglam,0.05)) %>% ungroup()

figname<-paste(pathFigures,"FigLambdaalpha.pdf",sep="")
pdf(width = wdthpdf, height = hghtpdf,
    file = figname)
p1<-ggplot(SummLam)+
  geom_line(aes(x=Alpha,y=medianLogLam,col=Group),size=ls)+ #median
  #geom_ribbon(aes(x=Alpha,ymin=q05,ymax=q95,col=Group,fill=Group),
  #            alpha=0.01,fill="grey50",size=ls/2,linetype=2)+ #quantiles
  geom_text(data=SummLam[SummLam$Alpha==1,],
            aes(y=medianLogLam,col=Group,label=Group),x=1.01,size=ts/3)+
  scale_color_manual(values=colsMSE[1:5])+
  labs(y="Group penalty")+
  #scale_y_log10()+
  theme_bw()+
  theme(axis.text.x=element_text(size=ts),
        axis.text.y=element_text(size=ts),
        axis.title.x=element_text(size=ts+2),
        axis.title.y=element_text(size=ts+2),
        legend.text=element_text(size=ts),
        legend.title=element_text(size=ts+2),
        strip.text=element_text(size=ts))#,
p1
dev.off()
p1

