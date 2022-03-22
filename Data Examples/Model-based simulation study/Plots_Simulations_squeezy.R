###Simulation study
#NOTE 1: in Rstudio; press ALT+o for overview of code section headings
#This file may be used to replicate the figures for Section 3.1.1-3.1.3

#Install and load libraries----
setwd("") #make sure working directory is set to where this file is
pathResults <- "" #if desired, separate folder can be specified to save results..
pathFigures <- "" #..and figures

#for storing and plotting results
if(!requireNamespace("dplyr")) install.packages("dplyr")
library(dplyr)
if(!requireNamespace("ggplot2")) install.packages("ggplot2")
library(ggplot2)
if(!requireNamespace("RColorBrewer")) install.packages("RColorBrewer")
library(RColorBrewer)
if(!requireNamespace("scales")) install.packages("scales")
library(scales)
if(!requireNamespace("tidyr")) install.packages("tidyr")
library(tidyr)
if(!requireNamespace("ggpubr")) install.packages("ggpubr")
library(ggpubr)
if(!requireNamespace("ggrepel")) install.packages("ggrepel")
library(ggrepel)
if(!requireNamespace("scales")) install.packages("scales")
library(scales)



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
colsAUC2 <- brewer.pal(10,"RdYlBu")
colsMSE <- seq_gradient_pal(low="grey50",high="white")((0:8)/9)
colsTime <- seq_gradient_pal(low="black",high="grey50")((0:8)/9)

#load data----
n<-150
p<-1200
#load data for plots performance MSE
alphaEN <- 0.3
df2 <- data.frame()
for(alphaEN in c(0.3,0.8)){
  for(setting in c("1","2","3")){
    #glmnet, fwelnet, iplasso (small grid)
    fname <- paste(pathResults,paste("simresEN",alphaEN,setting,sep="_"),".Rdata",sep="") 
    load(fname)
    df$Setting <- setting
    df$Alpha <- alphaEN
    df2 <- rbind(df2,df)
    
    #squeezy
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

#load data for plots time
df3 <- data.frame()
alphaEN <- 0.3
for(setting in as.character(4:11)){
  fname <- paste(pathResults,paste("simresEN",alphaEN,setting,sep="_"),".Rdata",sep="")
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


#load data for plot MSE for range alpha
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

#load data for plot group estimates
alphaEN <- 1
n<-150
p<-1200
dfGrps2 = data.frame()
for(setting in as.character(1:3)){
  fname <- paste(pathResults,paste("simresSqueezy",alphaEN,setting,sep="_"),".Rdata",sep="")
  load(fname)
  dfGrps$Sigma2[dfGrps$Method=="squeezy.m"] <- dfGrps$Sigma2[dfGrps$Method=="squeezy.m"]
  dfGrps$Lambda[dfGrps$Method=="squeezy.m"] <- dfGrps$Lambda[dfGrps$Method=="squeezy.m"]
  dfGrps$Tau2[dfGrps$Method=="squeezy.m"] <- dfGrps$Sigma2[dfGrps$Method=="squeezy.m"]/dfGrps$Lambda[dfGrps$Method=="squeezy.m"]
  
  dfGrps2 = rbind(dfGrps2,dfGrps)
  
  
  #empirical truth
  fname <- paste(pathResults,paste("simresSqueezy_GroupsEmpTrue",alphaEN,setting,sep="_"),".Rdata",sep="")
  load(fname)
  dfGrps2 = rbind(dfGrps2,dfGrps)
}
dfGrps2$Setting <- factor(dfGrps2$Setting,levels=c("1","2","3"),
                          labels=c("No groups","Weakly informative groups","Informative groups"))
dfGrps2$Method <- factor(dfGrps2$Method ,levels=unique(dfGrps2$Method)[c(2,1,4,3)],
                         labels=c("ecpcEN squeezy", "squeezy (multi)", "sample estimate", "Truth"))


#Plot: prediction MSE----
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

#Plot: prediction informative and reCV multiple alpha----
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

p2<-ggplot(df2[df2$Alpha==alphaEN & df2$Setting=="Informative groups",])+aes(x=MethodNum,y=MSE)+
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
        legend.position ="right")#,
p2

figname<-paste(pathFigures,"FigSimsPredAll.pdf",sep="")
pdf(width = wdthpdf*1.7, height = hghtpdf,
    file = figname)
p4<-ggarrange(p1,p2,common.legend=F,nrow=1, widths=c(1,1.3),legend="right")+
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



#Plot: group lambda estimates versus truth----

lb <- min(1/dfGrps2$Lambda[dfGrps2$Lambda<Inf])
dfGrps2$Lambdacutoff <- pmin(dfGrps2$Lambda,10^7)
figname<-paste(pathFigures,"FigSimsEstimatesLambda",alphaEN,".pdf",sep="")
pdf(width = wdthpdf*1.5, height = hghtpdf*1.1,
    file = figname)
p1<-ggplot(dfGrps2[dfGrps2$Method!="Truth",])+aes(x=factor(Group),y=Lambdacutoff)+
  geom_boxplot(aes(fill=Method))+
  facet_wrap(.~Setting,scales="free")+
  scale_y_log10()+
  #scale_fill_manual(values=colsMSE[1:3],name="Method")+
  labs(y=expression(lambda[g]), x="Group")+
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

figname<-paste(pathFigures,"FigSimsEstimatesTau",alphaEN,".pdf",sep="")
pdf(width = wdthpdf*1.5, height = hghtpdf*1.1,
    file = figname)
p1<-ggplot(dfGrps2[dfGrps2$Method!="Truth",])+aes(x=factor(Group),y=Tau2)+
  geom_boxplot(aes(fill=Method))+
  facet_wrap(.~Setting,scales="free")+
  #scale_y_log10()+
  #scale_fill_manual(values=colsMSE[1:3],name="Method")+
  labs(y=expression(tau[g]^2),x="Group")+
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
            aes(y=meanTime),linetype=2,size=ls,alpha=0.3)+
  #scale_color_manual(values=colsTime[1:9],name="Method")+
  scale_color_manual(values=colsAUC2[-5],name="Method")+
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
            aes(y=meanTime),linetype=2,size=ls,alpha=0.3)+
  #scale_color_manual(values=colsTime[1:9],name="Method")+
  scale_color_manual(values=colsAUC2[-5],name="Method")+
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
            aes(y=meanTime),linetype=2,size=ls,alpha=0.3)+
  #scale_color_manual(values=colsTime[1:9],name="Method")+
  scale_color_manual(values=colsAUC2[-5],name="Method")+
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
