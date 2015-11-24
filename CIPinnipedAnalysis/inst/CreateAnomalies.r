#if fdir="" it looks for data files in the CalcurData package directory of your R library.
#if fdir=NULL it looks in databases.txt in CalcurData package directory to get the database location
#if fdir is anything else it uses the value of fdir as the directory for database.  
#The scripts check for the value of fdir and if it exists the script will not change the value; otherwise it sets it to NULL
if(!exists("fdir"))fdir=NULL
require(CIPinnipedAnalysis)
if(!exists("lastyear"))lastyear=2015
if(!exists("locations"))locations=2:5
# use "" to use databases in Calcur installed package directory; use NULL to use default Databases directory J:/Master  or specify directory
#fdir=NULL
################################################################################
# Construct environmental values for models
#####################################################
# Sea level height data: http://tidesandcurrents.noaa.gov/sltrends/sltrends_station.shtml?stnid=9412110
#####################################################
data(SeaLevelHeight)
JunetoSeptSLH=with(SeaLevelHeight[SeaLevelHeight$Month%in% 6:9,],tapply(SeaLevelHeight,Year,mean,na.rm=T))
AprtoSeptSLH=with(SeaLevelHeight[SeaLevelHeight$Month%in% 4:9,],tapply(SeaLevelHeight,Year,mean,na.rm=T))
octtodec=split(SeaLevelHeight$SeaLevelHeight[SeaLevelHeight$Month%in%10:12],SeaLevelHeight$Year[SeaLevelHeight$Month%in%10:12])
jantofeb=split(SeaLevelHeight$SeaLevelHeight[SeaLevelHeight$Month%in%1:2],SeaLevelHeight$Year[SeaLevelHeight$Month%in%1:2])
OcttoFebSLH=NULL
for(i in 1:length(octtodec))
{
	if(i+1 > length(jantofeb))
	  OcttoFebSLH=c(OcttoFebSLH,mean(octtodec[[i]]))
    else
		OcttoFebSLH=c(OcttoFebSLH,mean(c(octtodec[[i]],jantofeb[[i+1]])))
}
numyears=lastyear-1975+1
JunetoSeptSLH=JunetoSeptSLH[1:numyears]
AprtoSeptSLH=AprtoSeptSLH[1:numyears]
OcttoFebSLH=OcttoFebSLH[1:numyears]
#####################################################
# SST data
#####################################################
# Create SST Anomalies
anomalies=create.SST.anomalies(1975:lastyear,fdir=fdir)

# Use Locations 2-5 (WSB,PtArg,PtSM,PtSL) for pup weight predictions
SSTAnomalies=t(apply(anomalies[,,locations],c(2,1),mean,na.rm=TRUE))
SSTAnomalies[is.nan(SSTAnomalies)]=NA
# Set maxyear, minyear and numyears,lastyear
maxyear= max(as.numeric(row.names(SSTAnomalies)))
minyear= min(as.numeric(row.names(SSTAnomalies)))
numyears=maxyear-minyear+1

AprtoSeptAnomalies=average_anomalies(SSTAnomalies,4,6)
JunetoSeptAnomalies=average_anomalies(SSTAnomalies,6,4)
JulytoJuneAnomalies=average_anomalies(SSTAnomalies,7,12)
JunetoFebAnomalies=average_anomalies(SSTAnomalies,6,9)
OcttoDecAnomalies=average_anomalies(SSTAnomalies,10,3)
OcttoFebAnomalies=average_anomalies(SSTAnomalies,10,5)
JantoMayAnomalies=average_anomalies(SSTAnomalies,1,5)
SSTAnomalyBySeason=rbind(data.frame(Year=minyear:maxyear,Season=rep("Spring",numyears),SSTAnomaly=as.vector(JantoMayAnomalies)),
		data.frame(Year=minyear:maxyear,Season=rep("Summer",numyears),SSTAnomaly=as.vector(JunetoSeptAnomalies)),
		data.frame(Year=minyear:maxyear,Season=rep("Fall",numyears),SSTAnomaly=as.vector(OcttoDecAnomalies)) )
SSTAnomalyBySeason$Season=factor(SSTAnomalyBySeason$Season,levels=c("Spring","Summer","Fall"))
SSTAnomalyBySeason=SSTAnomalyBySeason[order(SSTAnomalyBySeason$Year,SSTAnomalyBySeason$Season),]  # this is not used in this script yet


# Compute SST anomaly averages across month grouping to use with period specific growth rate
#JantoMayAnomalies=rowMeans(SSTAnomalies[,c("Jan","Feb","Mar","Apr","May")],na.rm=TRUE)[1:numyears]
#OcttoDec=SSTAnomalies[,c("Oct","Nov","Dec")][1:numyears,]
#JantoFeb=SSTAnomalies[2:nrow(SSTAnomalies),c("Jan","Feb")]
#JantoJune=SSTAnomalies[2:nrow(SSTAnomalies),c("Jan","Feb","Mar","Apr","May","June")]
#if(nrow(JantoFeb)<nrow(OcttoDec)) 
#{
#	JantoFeb=rbind(JantoFeb,c(NA,NA))
#}else
#{
#	if(nrow(JantoFeb)>nrow(OcttoDec)) JantoFeb=JantoFeb[-nrow(JantoFeb),]
#}
#
#
#OcttoFebAnomalies=as.matrix(cbind(OcttoDec,JantoFeb))
#OcttoFebAnomalies[is.nan(OcttoFebAnomalies)]=NA
#OcttoFebAnomalies=rowMeans(OcttoFebAnomalies,na.rm=TRUE)
#OcttoJuneAnomalies=as.matrix(cbind(OcttoDec,JantoJune))
#OcttoJuneAnomalies[is.nan(OcttoJuneAnomalies)]=NA
#OcttoJuneAnomalies=rowMeans(OcttoJuneAnomalies,na.rm=TRUE)
#
#JunetoSeptAnomalies=average_anomalies(SSTAnomalies,6,4)
#JulytoJuneAnomalies=average_anomalies(SSTAnomalies,7,12)
#JunetoFebAnomalies=average_anomalies(SSTAnomalies,6,9)
#OcttoDecAnomalies=average_anomalies(SSTAnomalies,10,3)
#JantoMayAnomalies=average_anomalies(SSTAnomalies,1,5)

#JunetoSeptAnomalies=SSTAnomalies[1:numyears,c("June","July","Aug","Sept")]
#JunetoSeptAnomalies[is.nan(JunetoSeptAnomalies)]=NA
#JunetoSeptAnomalies=rowMeans(JunetoSeptAnomalies,na.rm=TRUE)[1:numyears]
#AprtoSeptAnomalies=SSTAnomalies[1:numyears,c("Apr","May","June","July","Aug","Sept")]
#AprtoSeptAnomalies[is.nan(AprtoSeptAnomalies)]=NA
#AprtoSeptAnomalies=rowMeans(AprtoSeptAnomalies,na.rm=TRUE)[1:numyears]
#OcttoDecAnomalies=SSTAnomalies[1:numyears,c("Oct","Nov","Dec")]
#OcttoDecAnomalies[is.nan(OcttoDecAnomalies)]=NA
#OcttoDecAnomalies=rowMeans(OcttoDecAnomalies,na.rm=TRUE)[1:numyears]
#JunetoFebAnomalies=as.matrix(cbind(SSTAnomalies[1:numyears,c("June","July","Aug","Sept","Oct","Nov","Dec")],JantoFeb))
#JunetoFebAnomalies[is.nan(JunetoFebAnomalies)]=NA
#JunetoFebAnomalies=rowMeans(JunetoFebAnomalies,na.rm=TRUE)[1:numyears]
#JulytoJuneAnomalies=as.matrix(cbind(SSTAnomalies[1:numyears,c("July","Aug","Sept","Oct","Nov","Dec")],JantoJune))
#JulytoJuneAnomalies[is.nan(JulytoJuneAnomalies)]=NA
#JulytoJuneAnomalies=rowMeans(JulytoJuneAnomalies,na.rm=TRUE)[1:numyears]
#SSTAnomalyBySeason=rbind(data.frame(Year=minyear:lastyear,Season=rep("Spring",numyears),SSTAnomaly=as.vector(JantoMayAnomalies)),
#		data.frame(Year=minyear:lastyear,Season=rep("Summer",numyears),SSTAnomaly=as.vector(JunetoSeptAnomalies)),
#		data.frame(Year=minyear:lastyear,Season=rep("Fall",numyears),SSTAnomaly=as.vector(OcttoDecAnomalies)) )
#SSTAnomalyBySeason$Season=factor(SSTAnomalyBySeason$Season,levels=c("Spring","Summer","Fall"))
#SSTAnomalyBySeason=SSTAnomalyBySeason[order(SSTAnomalyBySeason$Year,SSTAnomalyBySeason$Season),]  # this is not used in this script yet

#####################################################
# Upwelling index
#####################################################
# Extract UWI data
UWI=getCalcurData("Environ","UWIAnomaly",dir=fdir)

UWI33=UWI[UWI$Location=="33N119W",]
UWI36=UWI[UWI$Location=="36N122W",]
rownames(UWI33)=UWI33$Year
rownames(UWI36)=UWI36$Year
UWI33=UWI33[,-(1:2)]
UWI36=UWI36[,-(1:2)]

UWImeansJunetoSept=rbind(average_anomalies(UWI33,6,4),average_anomalies(UWI36,6,4))
UWImeansAprtoSept=rbind(average_anomalies(UWI33,4,6),average_anomalies(UWI36,6,6))
UWImeansOcttoFeb=rbind(average_anomalies(UWI33,10,5),average_anomalies(UWI36,10,5))
UWImeansJunetoFeb=rbind(average_anomalies(UWI33,6,9),average_anomalies(UWI36,6,9))
UWImeansOcttoDec=rbind(average_anomalies(UWI33,10,3),average_anomalies(UWI36,10,3))

rownames(UWImeansJunetoSept)=c("33N119W","36N122W")
rownames(UWImeansAprtoSept)=c("33N119W","36N122W")
rownames(UWImeansOcttoFeb)=c("33N119W","36N122W")
rownames(UWImeansJunetoFeb)=c("33N119W","36N122W")
rownames(UWImeansOcttoDec)=c("33N119W","36N122W")


UWImeansJunetoSept=UWImeansJunetoSept[,as.numeric(colnames(UWImeansJunetoSept))<=lastyear]
UWImeansAprtoSept=UWImeansAprtoSept[,as.numeric(colnames(UWImeansAprtoSept))<=lastyear]
UWImeansOcttoFeb=UWImeansOcttoFeb[,as.numeric(colnames(UWImeansOcttoFeb))<=lastyear]
UWImeansJunetoFeb=UWImeansJunetoFeb[,as.numeric(colnames(UWImeansJunetoFeb))<=lastyear]

#UWI=UWI[order(UWI$Year,UWI$Month),]
#UWI=UWI[UWI$Year<=lastyear+1,]
#UWImeansJunetoSept=with(UWI[UWI$Month%in%6:9,], tapply(UWIAnomaly,list(Location,Year),mean,na.rm=TRUE))
#UWIJunetoSept=with(UWI[UWI$Month%in%6:9,], tapply(UWIAnomaly,list(Month,Year,Location),mean,na.rm=TRUE))
#UWImeansAprtoSept=with(UWI[UWI$Month%in%4:9,], tapply(UWIAnomaly,list(Location,Year),mean,na.rm=TRUE))
#UWIAprtoSept=with(UWI[UWI$Month%in%4:9,], tapply(UWIAnomaly,list(Month,Year,Location),mean,na.rm=TRUE))
#UWIOcttoDec=with(UWI[UWI$Month%in%10:12,], tapply(UWIAnomaly,list(Month,Year,Location),mean,na.rm=TRUE))
#UWIJantoFeb=with(UWI[UWI$Month%in%1:2,], tapply(UWIAnomaly,list(Month,Year,Location),mean,na.rm=TRUE))

#minyr=min(dim(UWIJunetoSept)[2],dim(UWIOcttoDec)[2],dim(UWIJantoFeb)[2]-1)

#UWIAprtoSept=UWIAprtoSept[,1:minyr,]
#UWIJunetoSept=UWIJunetoSept[,1:minyr,]
#UWIOcttoDec=UWIOcttoDec[,1:minyr,]

#UWImeansOcttoFeb=NULL
#for(i in 1:2)
#{
#	if(dim(UWIOcttoDec)[2]>=dim(UWIJantoFeb)[2])
#		UWImeansOcttoFeb=rbind(UWImeansOcttoFeb,colMeans(rbind(UWIOcttoDec[,-dim(UWIOcttoDec)[2],i],UWIJantoFeb[,-1,i]),na.rm=TRUE))
#	else
#		UWImeansOcttoFeb=rbind(UWImeansOcttoFeb,colMeans(rbind(UWIOcttoDec[,,i],UWIJantoFeb[,-1,i]),na.rm=TRUE))
#}
#UWImeansJunetoFeb=NULL
#for(i in 1:2)
#{
#	if(dim(UWIJunetoSept)[2]<=dim(UWIOcttoDec)[2])
#		UWImeansJunetoFeb=rbind(UWImeansJunetoFeb,colMeans(rbind(UWIJunetoSept[,,i],UWIOcttoDec[,,i],UWIJantoFeb[,-1,i]),na.rm=TRUE))
#    else
#		if(dim(UWIOcttoDec)[2]>dim(UWIJantoFeb)[2])
#			UWImeansJunetoFeb=rbind(UWImeansJunetoFeb,colMeans(rbind(UWIJunetoSept[,-dim(UWIJunetoSept)[2],i],UWIOcttoDec[,-dim(UWIOcttoDec)[2],i],UWIJantoFeb[,-1,i]),na.rm=TRUE))
#        else
#			UWImeansJunetoFeb=rbind(UWImeansJunetoFeb,colMeans(rbind(UWIJunetoSept[,-dim(UWIJunetoSept)[2],i],UWIOcttoDec[,,i],UWIJantoFeb[,-1,i]),na.rm=TRUE))
#}
#UWImeansJunetoSept=UWImeansJunetoSept[,as.numeric(colnames(UWImeansJunetoSept))<=lastyear]
#UWImeansAprtoSept=UWImeansAprtoSept[,as.numeric(colnames(UWImeansAprtoSept))<=lastyear]
#UWImeansOcttoFeb=UWImeansOcttoFeb[,as.numeric(colnames(UWImeansOcttoFeb))<=lastyear]
#UWImeansJunetoFeb=UWImeansJunetoFeb[,as.numeric(colnames(UWImeansJunetoFeb))<=lastyear]

####################################################
# Multivariate ENSO Index
####################################################
# Extract MEI data
MEI=getCalcurData("Environ","MEI",dir=fdir)
MEI=MEI[order(MEI$Year,MEI$Month),]
# Compute correlations between MEI and SST to find the best lag to use for MEI
SSTAnomalies.db=data.frame(SSTAnomaly=as.vector(t(SSTAnomalies[-(1:2),])))
MEIcor=vector("numeric",8)
for(lag in 0:7){
	MEIcor[lag+1]=cor(MEI$MEI[1:(length(MEI$MEI)-lag)],SSTAnomalies.db$SSTAnomaly[(lag+1):length(MEI$MEI)],use="complete.obs")
	cat("\nlag = ",lag,"cor = ",MEIcor[lag+1])
}
lag=which(MEIcor==max(MEIcor))-1
# lag is 2 months - calculation for OcttoFeb assumes 2 month lag
average.MEI=function(x,months)return(tapply(x$MEI[x$Month%in%months],x$Year[x$Month%in%months],mean))
LaggedMEIJunetoSept=average.MEI(MEI,(6:9-lag))
LaggedMEIJunetoSept=LaggedMEIJunetoSept[as.numeric(names(LaggedMEIJunetoSept))<=lastyear]
LaggedMEIAprtoSept=average.MEI(MEI,(4:9-lag))
LaggedMEIAprtoSept=LaggedMEIAprtoSept[as.numeric(names(LaggedMEIAprtoSept))<=lastyear]

# next assume 2 month lag to avoid Dec/Jan break
LaggedMEIOcttoFeb=average.MEI(MEI,8:12) 
LaggedMEIOcttoFeb=LaggedMEIOcttoFeb[as.numeric(names(LaggedMEIOcttoFeb))<=lastyear]

LaggedMEIJunetoFeb=average.MEI(MEI,4:12)
LaggedMEIJunetoFeb=LaggedMEIJunetoFeb[as.numeric(names(LaggedMEIJunetoFeb))<=lastyear]


