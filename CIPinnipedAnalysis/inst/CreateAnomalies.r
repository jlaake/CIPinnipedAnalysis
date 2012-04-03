anomalies=create.SST.anomalies(c(1994:1996,1998:2008))
CentralSSTAnomalies=t(apply(anomalies[,,3:7],c(2,1),mean,na.rm=TRUE))
SouthSSTAnomalies=t(apply(anomalies[,,1:2],c(2,1),mean,na.rm=TRUE))
NorthSSTAnomalies=anomalies[,,8]
SSTAnomalies=t(apply(anomalies[,,1:5],c(2,1),mean,na.rm=TRUE))
maxyear= max(as.numeric(row.names(SSTAnomalies)))
minyear= min(as.numeric(row.names(SSTAnomalies)))
numyears=maxyear-minyear+1
JantoMayAnomalies=rowMeans(SSTAnomalies[,c("Jan","Feb","Mar","Apr","May")])
JunetoSeptAnomalies=rowMeans(SSTAnomalies[,c("June","July","Aug","Sept")])
OcttoDecAnomalies=rowMeans(SSTAnomalies[,c("Oct","Nov","Dec")])
x=rbind(data.frame(Year=minyear:maxyear,Season=rep("Spring",numyears),SSTAnomaly=JantoMayAnomalies),
		data.frame(Year=minyear:maxyear,Season=rep("Summer",numyears),SSTAnomaly=JunetoSeptAnomalies),
		data.frame(Year=minyear:maxyear,Season=rep("Fall",numyears),SSTAnomaly=OcttoDecAnomalies) )
x$Season=factor(x$Season,levels=c("Spring","Summer","Fall"))
x=x[order(x$Year,x$Season),]
pdf("SSTAnomaly.pdf",pointsize=10)
par(mfrow=c(3,1))
plot(minyear:maxyear,x$SSTAnomaly[x$Season=="Spring"],ylab="SST Anomaly",main="Winter/Spring (Jan-May)",xlab="")
lines(minyear:maxyear,x$SSTAnomaly[x$Season=="Spring"])
abline(h=0)
plot(minyear:maxyear,x$SSTAnomaly[x$Season=="Summer"],ylab="SST Anomaly",main="Summer (June-Sept)",xlab="")
lines(minyear:maxyear,x$SSTAnomaly[x$Season=="Summer"])
abline(h=0)
plot(minyear:maxyear,x$SSTAnomaly[x$Season=="Fall"],ylab="SST Anomaly",main="Fall (Oct-Dec)",xlab="")
lines(minyear:maxyear,x$SSTAnomaly[x$Season=="Fall"])
abline(h=0)
dev.off()

pdf("MultivariateENSOIndex.pdf",pointsize=10)
MEI=getCalcurData("Environ","MEI")
minyear=min(MEI$Year)
maxyear=max(MEI$Year)
numyears=maxyear-minyear+1
plot(MEI$MEI,type="l",lwd=2,xaxt="n",ylab="MEI",xlab="Year")
axis(1,at=12*(0:(numyears-1))+1,labels=as.character(minyear:maxyear))
abline(h=0)
abline(h=1)
abline(h=-1)
dev.off()

pdf("UWIAnomaly.pdf",pointsize=10)
par(mfrow=c(2,1))
UWI=getCalcurData("Environ","UWIAnomaly")
UWI=UWI[order(UWI$Year,UWI$Month),]
minyear=min(UWI$Year)
maxyear=max(UWI$Year)
numyears=maxyear-minyear+1
plot(UWI$UWI[UWI$Location=="36N122W"],type="l",lwd=2,xaxt="n",ylab="UWI",xlab="Year",main="36N122W")
axis(1,at=12*(0:(numyears-1))+1,labels=as.character(minyear:maxyear))
abline(h=0)
plot(UWI$UWI[UWI$Location=="33N119W"],type="l",lwd=2,xaxt="n",ylab="UWI",xlab="Year",main="33N119W")
axis(1,at=12*(0:(numyears-1))+1,labels=as.character(minyear:maxyear))
abline(h=0)
dev.off()


