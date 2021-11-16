# OrderYourMind_fMRI
# =============================================================
#Removing evrything from workspace
graphics.off()
rm(list = ls(all = TRUE))

# Load libraries
library(data.table)
library(foreign)
library(dplyr)
library(tidyr)
library(gridExtra)
library(ggplot2)
library(RColorBrewer)
library(corrplot)
library(gee)
library(broom)
library(effects)
library(GGally)
library(psych)
library(stargazer)
library(afex)
library(lme4)

#Setting up directory (e.g., directories for AK)
fsep<-.Platform$file.sep;
Dir_R<-path.expand("~/Dropbox/w_ONGOINGRFILES/w_OTHERS")
Dir_BDATA<-path.expand("~/Dropbox/sharedprojects/P.SEQ_fMRI/DataAnalysis/MATLAB/BEH_S2MRI/w_ALLGRAND")

# Source Files
setwd(Dir_R)
source('basic_theme.R')
source('basic_libS.R')
# source('anovakun_480.txt')
#======================================================================================================
# # Load data
#======================================================================================================

# READING FILES
setwd(Dir_BDATA)#Change directory
ds<-fread('OrderYourMind_fMRI_BEHP.txt', header=TRUE, fill=TRUE)

#======================================================================================================
# # Edit & Select data 
#======================================================================================================

# Add variables
ds[is.nanM(ds)]<-NA;
ds$RTlog<-log(ds$RT);
ds[,ACPOS:=if_else(BSPOS<=3,1,2)]
ds[,SESSION:=if_else(SUBID>=200,"fMRI","Prep")]
ds[SUBID>=200,SUBID:=SUBID-100]

# Recode
ds$WCPOS_L<-dplyr::recode(ds$WCPOS,`1`="1st",`2`="2nd",`3`="3rd")
ds$ACPOS_L<-dplyr::recode(ds$ACPOS,`1`="1st",`2`="2nd")

# Factorize
cols<-c("WCPOS_L","ACPOS_L","LAG2","BSPOS")
ds[,(cols):=lapply(.SD,factor),.SDcols=cols];

#Selection Rule!!
ds<-ds[TRIAL>1]
ds_p<-ds[PRACTICE==0]#Practice
ds<-ds[PRACTICE==1]#Experiment (include both prep and fMRI)
# ds<-subset(ds,SUBID!=104)

# RT cutoff
# cutRT<-quantile(ds$RT[ds$ACC==1&ds$LACC==1],0.995,na.rm=T);
# print(paste0("0.5% Response cutoff: ", round(cutRT,2)));
# ds<-ds[ds$RT>10 & ds$RT<cutRT,];
ds_f<-subset(ds,ACC==1 & LACC==1);

# Checking
rt<-ds_f %>% group_by(SESSION,SUBID) %>% summarise(DV=mean(RT,na.rm=TRUE));as.data.frame(rt)
acc<-ds %>% group_by(SESSION,SUBID) %>% summarise(DV=mean(ACC,na.rm=TRUE));as.data.frame(acc)
# s_s<-rt %>% tidyr::spread(DV);plot(s_s$Session1,s_s$Session2)

#======================================================================================================
# # Aggregation
#======================================================================================================

# RT:BSPOS
bspos_RT<-ds_f%>%group_by(SUBID,SESSION,BSPOS)%>%summarise(DV=mean(RT,na.rm=TRUE))%>%
  normWS("SUBID","DV")%>%group_by(SESSION,BSPOS)%>%summarize(DVm=mean(DV),se=sd(DV_n)/sqrt(n()),wseb=se*1.96,n=n());print(bspos_RT)

# ACC:BSPOS
bspos_ACC<-ds%>%group_by(SUBID,SESSION,BSPOS)%>%summarise(DV=(1-mean(ACC,na.rm=TRUE))*100)%>%
  normWS("SUBID","DV")%>%group_by(SESSION,BSPOS)%>%summarize(DVm=mean(DV),se=sd(DV_n)/sqrt(n()),wseb=se*1.96,n=n());print(bspos_ACC)

# RT: Learning
learn_RT<-ds_f%>%group_by(SUBID,BLOCK)%>%summarise(DV=mean(RT,na.rm=TRUE))%>%
  normWS("SUBID","DV")%>%group_by(BLOCK)%>%summarize(DVm=mean(DV),se=sd(DV_n)/sqrt(n()),wseb=se*1.96,n=n());print(learn_RT)

# ACC: Learning
learn_ACC<-ds%>%group_by(SUBID,BLOCK)%>%summarise(DV=mean(ACC,na.rm=TRUE))%>%
  normWS("SUBID","DV")%>%group_by(BLOCK)%>%summarize(DVm=mean(DV),se=sd(DV_n)/sqrt(n()),wseb=se*1.96,n=n());print(ds_f)

# Individuals
bspos_RTi<-ds_f%>%group_by(SUBID,BSPOS)%>%summarise(DV=mean(RT,na.rm=TRUE))
bspos_ACCi<-ds%>%group_by(SUBID,BSPOS)%>%summarise(DV=mean(ACC,na.rm=TRUE))


#======================================================================================================
# # Individual difference!
#======================================================================================================


# #======================================================================================================
# #Statistical tests 
# #======================================================================================================

# # Switch effect
# m_rt=lmer(RTlog~1+SWITCH*SESSION+(1+SWITCH*SESSION|SUBID),data=ds_f);summary(m_rt)
# m_acc=glmer(ACCs~1+SWITCH*SESSION+(1+SWITCH*SESSION|SUBID),family=binomial,data=ds);summary(m_acc)
# 
# # Switch effect (aggregated levels)
# aggRT<-ds_f %>% group_by(SUBID,SWITCH) %>% summarise(DV=mean(RT,na.rm=TRUE))
# aggACC<-ds %>% group_by(SUBID,SWITCH) %>% summarise(DV=mean(ACC,na.rm=TRUE))
# r_rt_bi<-effS(ezANOVA(aggRT, DV,SUBID, within = c(SWITCH), return_aov=T,detailed=T));print(r_rt_bi$ANOVA)
# r_acc_bi<-effS(ezANOVA(aggACC, DV,SUBID, within = c(SWITCH),return_aov=T,detailed=T));print(r_acc_bi$ANOVA)
# 

# model.matrix(m_it)
#======================================================================================================
#Process3-Plotting-ggplot2
#======================================================================================================
theme_set(theme_bw(base_size = 20))#32/28
CSCALE_PURD = rev(brewer.pal(9,"PuRd"));
CSCALE_BLUE = rev(brewer.pal(9,"Blues"));
CSCALE_PiYG = rev(brewer.pal(11,"PiYG"));
CSCALE_RdBu = rev(brewer.pal(11,"RdBu"));
CSCALE_PAIRED = rev(brewer.pal(12,"Paired"));
CSCALE_YlGnBu = rev(brewer.pal(9,"YlGnBu"));
CSCALE_BrBG = rev(brewer.pal(9,"BrBG"));
CSCALE_Set1 = (brewer.pal(9,"Set1"));
CSCALE_Pair = (brewer.pal(9,"Paired"));
CSCALE_Grey = (brewer.pal(9,"Greys"));
CSCALE_Accent = (brewer.pal(8,"Accent"));
CSETS<-CSCALE_Accent[5:8]
CSETS<-c("red",'black')

# Setting
actDS<-subset(bspos_RT,SESSION=="fMRI")

quartz(width=5.5,height=4.70)
ggplot(actDS, aes(x=BSPOS,y=DVm,group=1,color="black")) + 
  geom_line(size=1.5) +
  geom_point(size=5) +
  #facet_wrap(~SESSION)+
  #coord_cartesian(ylim=c(325,425))+
  #coord_cartesian(ylim=c(-1,8))+
  scale_color_manual(values=CSETS) +
  geom_errorbar(aes(ymin=DVm-wseb, ymax=DVm+wseb),width=.15,size=1)+#,color="black"
  ylab("Response Time(ms)")+
  #ylab("Error(p)")+
  xlab("Sequence Positions")+
  #Aesthetics!-------------------------
  theme(panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.margin=unit(rep(1,1,4),"line"),
      axis.line=element_line(size=1),
      legend.key = element_blank(),
      #legend.key.size=unit(1,"lines"),
      legend.position = "none",
      legend.title=element_blank(),
      legend.text =element_text(family="Helvetica",size=8),
      axis.ticks.length = unit(0.5, "lines"),
      axis.title  = element_text(family="Helvetica", face="bold",vjust=1.8),
      axis.text   = element_text(size=16),
      #axis.text.y = element_blank(),
      plot.title = element_blank(),
      strip.text = element_blank(),
      strip.background=element_blank())
