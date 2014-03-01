require(CIPinnipedAnalysis)
sdir=system.file(package="CIPinnipedAnalysis")
fdir=""
if(!exists("nboot"))nboot=100
lastyear=2013
# use "" to use databases in package directory; use NULL to use default Databases directory J:/Master  or specify directory
#dir=NULL

source(file.path(sdir,"CreateAnomalies.r"))
source(file.path(sdir,"ZC Weight Adjustment Model.r"))
source(file.path(sdir,"ZC Weight Environment Models.r"))
source(file.path(sdir,"ZC Growth.r"))

################################################################################
#  Brand vs control at brand eval
################################################################################
zcweights.be=get.zc.weights(fdir=fdir)
#
#  Use SMI from 1 Dec to end of Feb; use brand and unmarked only
#
zcweights.be=zcweights.be[zcweights.be$days>60&zcweights.be$days<150&zcweights.be$sitecode=="SMI"&zcweights.be$tag.type!="Tag",]
zcweights.be$tag.type=factor(zcweights.be$tag.type)
# lme model of tag.type effect
mod=lme(weight~tag.type*sex,random=~1|cohort,data=zcweights.be,method="ML")
summary(mod)
mod0=lme(weight~sex,random=~1|cohort,data=zcweights.be,method="ML")
summary(mod0)
mod1=lme(weight~tag.type+sex,random=~1|cohort,data=zcweights.be,method="ML")
summary(mod1)
anova(mod0,mod1,mod)
# kruskal test by sex
kruskal.test(weight~tag.type*cohort,data=zcweights.be[zcweights.be$sex=="F",])
kruskal.test(weight~tag.type*cohort,data=zcweights.be[zcweights.be$sex=="M",])
# paired t-test
be.fem=with(zcweights.be,tapply(weight,list(cohort,tag.type,sex),mean))[,,"F"]
be.m=with(zcweights.be,tapply(weight,list(cohort,tag.type,sex),mean))[,,"M"]
t.test(x=c(be.fem[,1],be.m[,1]),y=c(be.fem[,2],be.m[,2]),paired=TRUE)


zcweights.all=zcweights[zcweights$days>-30&zcweights$days<80&zcweights$cohort>=1997&
				 zcweights$sitecode=="SMI"&zcweights$type=="Initial",]
# Distribution of timing of weight samples
batches=apply(table(zcweights.all$cohort,zcweights.all$days),2,function(x) length(x[x>0]))
barplot(batches,xlab="Days from 1 Oct",ylab="Frequency",space=0)
abline(v=which(as.numeric(names(batches))==0),lwd=2)


par(mfrow=c(2,1))
batches=with(zcweights.all[zcweights.all$sex=="M",],aggregate(weight,list(cohort,days),mean))
names(batches)=c("Year","Days","Weight")
plot(batches$Year,batches$Weight,xlab="Year",ylab="Average male weight in batch")
batches=with(zcweights.all[zcweights.all$sex=="F",],aggregate(weight,list(cohort,days),mean))
names(batches)=c("Year","Days","Weight")
plot(batches$Year,batches$Weight,xlab="Year",ylab="Average female weight in batch")

