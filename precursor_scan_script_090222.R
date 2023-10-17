if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("mzR",force=TRUE)
BiocManager::install("msdata",force=TRUE)
library(mzR)
library(Rcpp)
library(msdata)
library(tidyverse)
library(stringr)
setwd("C:/Users/atran/OneDrive - University of Maryland Baltimore/JWJones_Lab Share/Anh Tran/R/tandem_MS_analysis")

mzxml <- system.file("_F9_BTLE_BH_nESI_1mMLiOH_MS2_051822_001.mzML",
                     package = "msdata")
aa <- openMSfile("_F9_BTLE_BH_nESI_1mMLiOH_MS2_051822_001.mzML")
runInfo(aa)
instrumentInfo(aa)
peaks(aa,1) ->scan1
write.csv(scan1,file="test1.csv")
scan1[scan1<265 & scan1>264]
  
plot(pl1[,1],pl1[,2],type="h",lwd=1)
chromatogram(aa)
header(aa,1005)->headerscan1
header(centroidspectrum,1)->headerscan2
headerscan2[3]->msLevel
headerscan2[15]->precursormz

openMSfile("_F9_LTLE_BH_nESI_1mMLiOH_MS2_051822_001.MZML") ->centroidspectrum
runInfo(centroidspectrum)

as_tibble(peaks(centroidspectrum,1))->scan1centroid
header(centroidspectrum)->headercentroid
as_tibble(peaks(centroidspectrum,9))->scan9centroid
scan1centroid$mz[scan1centroid>=264.26&scan1centroid<=264.28]
length(scan1centroid$mz[scan1centroid$mz>=264.26&scan1centroid$mz<=264.28])

filter(scan1centroid,mz<400 & mz>320 &
         abs(360-mz)<=min(abs(360-(scan1centroid[scan1centroid<400&scan1centroid>320]))))->tibbleA
#filter to get the closest row to mz value (360) given a mz interval (320-400)

tibble(1,2,3)->tibbleC
list(1,2,3,4)->listA

for(e in seq_along(1:4)){
paste0("_LTLE_BH_5x_ItMS2_00",e,"_053122_001.mzML")->filename
  
openMSfile(filename) ->centroidspectrum
tibble(mz=1,intensity=1)->tibbleA
c(1)->v1
c(1)->v2
for(i in seq_along(1:runInfo(centroidspectrum)[[1]])){
  header(centroidspectrum,i)[[3]]->v1[i]
  header(centroidspectrum,i)[[15]]->v2[i]
  184.07->product_ion
  0.02 ->mz_error
  as_tibble(peaks(centroidspectrum,i))->df1
  if (length(df1$mz[df1$mz>=(product_ion-mz_error)&df1$mz<=(product_ion+mz_error)]) == 0 | is.na(v2[i])){
    NA->tibbleA[i,]
  }else{
    filter(df1, mz>(product_ion-mz_error) & mz<(product_ion+mz_error) &
             abs(product_ion-mz)<=min(abs(product_ion-(df1[df1>(product_ion-mz_error)&df1<(product_ion+mz_error)]))))->tibble1
    add_row(tibbleA,tibble1)->tibbleA 
  }
}
tibble(v1,v2)%>% mutate(tibbleA)->tibbleB
#tibble of 1) ms level 2) precursor m/z 3)
tibbleB%>%group_by(v2)%>%summarise(sum_int=sum(intensity,na.rm=TRUE))%>%mutate(tibbleB%>%count(v2))%>%
filter(sum_int!=0)%>%filter(floor(v2)%%2==1)->tibbleC
#average intensity, adds number of scans(n), filters out mz without desired fragment, filters out non-even masses 
tibbleC->listA[[e]]
}

rbind(listA[[1]],listA[[2]],listA[[3]],listA[[4]])->tibbleE1

write.csv(tibbleE1,file=paste("HILIC_SM.csv"))



#function for precursor filter based on neutral loss
#filename- "file_name.mzML", nloss_mz: desired neutral loss (NOT fixed production ion)
#mz_error: range to search product ion within
183.0659->nloss_mz
0.05 -> mz_error
tibble(1,2,3)->tibbleC
list(1,2,3,4)->listB
for(e in seq_along(1:4)){
  paste0("_LTLE_BH_5x_ItMS2_00",e,"_053122_001.mzML")->filename
  openMSfile(filename) ->centroidspectrum
tibble(mz=1,intensity=1)->tibbleA
c(1)->v1
c(1)->v2
for(i in seq_along(1:runInfo(centroidspectrum)[[1]])){
  header(centroidspectrum,i)[[3]]->v1[i]
  header(centroidspectrum,i)[[15]]->v2[i]
  v2[i]-nloss_mz->m_minus_nloss
  as_tibble(peaks(centroidspectrum,i))->df1
  if ((length(df1$mz[df1$mz>=(m_minus_nloss-mz_error)&df1$mz<=(m_minus_nloss+mz_error)])==0)|is.na(m_minus_nloss)){
    NA->tibbleA[i,]
  }else{
    filter(df1,mz>(m_minus_nloss-mz_error) & mz<(m_minus_nloss+mz_error) &
             abs(m_minus_nloss-mz)<=min(abs(m_minus_nloss-(df1[df1>(m_minus_nloss-mz_error)&df1<(m_minus_nloss+mz_error)]))))->tibble1
    add_row(tibbleA,tibble1)->tibbleA 
  }
}
tibble(v1,v2)%>% mutate(tibbleA)->tibbleB
#tibble of 1) ms level 2) precursor m/z 3)

tibbleB2 <- na.omit(tibbleB)
tibbleB2%>%group_by(v2)%>%summarise(sum_int=sum(intensity,na.rm=TRUE))%>%mutate(tibbleB2%>%count(v2))%>%
  filter(sum_int!=0)%>%filter(floor(v2)%%2==0)->tibbleD

tibbleB%>%group_by(v2)%>%summarise(sum_int=sum(intensity,na.rm=TRUE))%>%mutate(tibbleB%>%count(v2))%>%
  filter(sum_int!=0)%>%filter(floor(v2)%%2==1)->tibbleC
#aumintensity, adds number of scans(n), filters out mz without desired fragment, filters out non-even masses 
tibbleC%>%mutate(n2=tibbleD[[3]])->tibbleC
tibbleC->listB[[e]]
}
rbind(listB[[1]],listB[[2]],listB[[3]],listB[[4]])->tibbleE2

precursor_filter_nloss("_F4_BTLE_BH_nESI_1mMLiOH_MS2_051822_001.mzML",416.29,0.02)

#Cer Acyl chain rules:
# 264 product,( precursor -264 -50 )/14 for saturated n-acyl chains- H+ adduct
#SM 184 product
#GlcCer is 264 and 180 NL
tibble(1,"264","180NL")->tibbleF
tibble(1,2,3)->filter1
rbind(tibbleE1[1],tibbleE2[1])%>% unique()->tibbleF

for(i in seq_along(1:length(tibbleF[[1]]))){
  if(tibbleF[[1]][i]%in% tibbleE1[[1]])
  {tibbleE1%>%filter(tibbleE1[[1]]==tibbleF[[1]][i]) ->filter1
  filter1[[2]]->tibbleF[i,2]
  }else{NA->tibbleF[i,2]}
}
for(i in seq_along(1:length(tibbleF[[1]]))){
  if(tibbleF[[1]][i]%in% tibbleE2[[1]])
  {tibbleE2%>%filter(tibbleE2[[1]]==tibbleF[[1]][i]) ->filter1
    filter1[[2]]->tibbleF[i,3]
  }else{NA->tibbleF[i,3]}
}
write.csv(tibbleF,file="tibbleCerC18.csv")->tibbleE2







rep(1,85)
tibble(tibbleE2[1],tibbleE2[2],col=rep(c(1,2,1,2,1),17))->tibbleplot
tibbleE2[2]

tibbleplot2%>%ggplot()+geom_linerange(aes(x=x,y=y,ymax=y,ymin=0,color=factor(z)))+
scale_color_manual(values=c("0"= "Green","1"="Red","2"="Black","3"="Blue","4"="Orange"))+
scale_y_continuous(breaks = seq(0,6000, by = 1000),limits=c(0,10000),expand = c(0, 0))+
scale_x_continuous(breaks = seq(0,2000, by = 50),limits=c(600,800))+
theme(panel.background=element_rect(fill="white"),
      panel.grid.major=element_line(colour="white"),
      panel.grid.minor=element_blank(),
      axis.line=element_line(colour="black"),
      axis.title.x=element_text(face="italic"))+
labs(color="Legend Title")+
xlab("m/z")+
ylab("Intensity")

?labs
"x"->colnames(tibbleplot[1])
col.na

tibble(x=tibbleplot[[1]],y=tibbleplot[[2]],z=tibbleplot[[3]])->tibbleplot2

tibbleplot2%>%ggplot()+geom_linerange(aes(x=x,y=y,ymax=y,ymin=0,color=factor(z)))+
  scale_color_manual(values=c("0"= "Green","1"="Red","2"="Black","3"="Blue","4"="Orange"))+
  scale_y_continuous(breaks = seq(0,6000, by = 1000),limits=c(0,6000),expand = c(0, 0))+
  scale_x_continuous(breaks = seq(0,2000, by = 200),limits=c(600,2000))+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major=element_line(colour="white"),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),
        axis.title.x=element_text(face="italic"))+
  labs(color="Legend Title")+
  xlab("m/z")+
  ylab("Intensity")

#m/z external calibration 
read_csv("raw_LTLEBH_F2_0mMLiOH_063022.csv")  ->df1
lm(formula= df1$`Center X` ~ df1$Height ,data=df1) ->df1fit
coef(df1fit)



df1%>% ggplot(aes(x=df1$`Center X`,y=df1$Height))+geom_point()+
  geom_abline(aes(intercept=coef(df1fit)[1],slope=coef(df1fit)[2]))

