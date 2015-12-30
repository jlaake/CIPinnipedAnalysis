#if fdir="" it looks for data files in the CalcurData package directory of your R library.
#if fdir=NULL it looks in databases.txt in CalcurData package directory to get the database location
#if fdir is anything else it uses the value of fdir as the directory for database.  
#The scripts check for the value of fdir and if it exists the script will not change the value; otherwise it sets it to NULL
if(!exists("fdir"))fdir=NULL
require(CIPinnipedAnalysis)
sdir=system.file(package="CIPinnipedAnalysis")
if(!exists("nboot"))nboot=100
source(file.path(sdir,"CreateAnomalies.r"))
cat("\nWeight adjust model\n")
source(file.path(sdir,"ZC_Weight_Adjustment_Model.r"))
cat("\nWeight environment model\n")
source(file.path(sdir,"ZC_Weight_Environment_Model.r"))
cat("\nWeight environment/fish model\n")
source(file.path(sdir,"ZC_Weight_FishBiomass_Model.r"))
cat("\nWeight environment/fish/diet model\n")
source(file.path(sdir,"ZC_Weight_Diet_Model.r"))
cat("\nLongitudinal growth model\n")
source(file.path(sdir,"ZC_Growth.r"))

#
# store ZcWeight.df in CIPinnipedCensusQuery
#
store_weights(ZCWeight.df,fdir=fdir)

################################################################################
#  Brand vs control at brand eval
################################################################################
zcweights.be=get.zc.weights(fdir=fdir)
zcweights.be=zcweights.be[zcweights.be$cohort<=lastyear,]

#
#  Use SMI from 1 Dec to end of Feb; use brand and unmarked only
#
zcweights.be=zcweights.be[zcweights.be$days>60&zcweights.be$days<150&zcweights.be$sitecode=="SMI"&zcweights.be$tag.type!="Tag",]
zcweights.be$tag.type=factor(zcweights.be$tag.type)
zcweights.be$cohort=factor(zcweights.be$cohort)
# lme model of tag.type effect
mod=lme(weight~tag.type*sex,random=~1|cohort,data=zcweights.be,method="ML")
summary(mod)
mod0=lme(weight~sex,random=~1|cohort,data=zcweights.be,method="ML")
summary(mod0)
mod1=lme(weight~tag.type+sex,random=~1|cohort,data=zcweights.be,method="ML")
summary(mod1)
anova(mod0,mod1,mod)
# wilcox test by sex
female_vecs=tapply(zcweights.be$weight[zcweights.be$sex=="F"],list(zcweights.be$cohort[zcweights.be$sex=="F"],zcweights.be$tag.type[zcweights.be$sex=="F"]),mean)
wilcox.test(female_vecs[,1],female_vecs[,2],paired=TRUE)
male_vecs=tapply(zcweights.be$weight[zcweights.be$sex=="M"],list(zcweights.be$cohort[zcweights.be$sex=="M"],zcweights.be$tag.type[zcweights.be$sex=="M"]),mean)
wilcox.test(male_vecs[,1],male_vecs[,2],paired=TRUE)
# combined wilcox
all_vecs=rbind(female_vecs,male_vecs)
wilcox.test(all_vecs[,1],all_vecs[,2],paired=TRUE)

# all paired t-test
be.fem=with(zcweights.be,tapply(weight,list(cohort,tag.type,sex),mean))[,,"F"]
be.m=with(zcweights.be,tapply(weight,list(cohort,tag.type,sex),mean))[,,"M"]
t.test(x=c(be.fem[,1],be.m[,1]),y=c(be.fem[,2],be.m[,2]),paired=TRUE)


zcweights.all=zcweights[zcweights$days>-30&zcweights$days<80&zcweights$cohort>=1997&
				 zcweights$sitecode=="SMI"&zcweights$type=="Initial",]
# Distribution of timing of weight samples
batches=apply(table(zcweights.all$cohort,zcweights.all$days),2,function(x) length(x[x>0]))
barplot(batches,xlab="Days from 1 Oct",ylab="Frequency",space=0)
abline(v=which(as.numeric(names(batches))==0),lwd=2)


jpeg("ZCWeightsByBatch.jpg",height=600,width=600,quality=100,pointsize=12)
par(mfrow=c(2,1))
batches=with(zcweights.all[zcweights.all$sex=="M",],aggregate(weight,list(cohort,days),mean))
names(batches)=c("Year","Days","Weight")
plot(batches$Year,batches$Weight,xlab="Year",ylab="Average male weight in batch")
batches=with(zcweights.all[zcweights.all$sex=="F",],aggregate(weight,list(cohort,days),mean))
names(batches)=c("Year","Days","Weight")
plot(batches$Year,batches$Weight,xlab="Year",ylab="Average female weight in batch")
dev.off()
