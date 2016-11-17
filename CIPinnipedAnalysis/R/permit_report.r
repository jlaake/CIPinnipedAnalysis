#' Create Zc Permit report 
#' 
#' Produces a table with tally of sampled animals for permit 
#' 
#' @param year numeric year for permit work
#' @param dir directory for location of files; dir=NULL uses default value in databases.txt in CalcurData
#' @export 
#' @author Jeff Laake
#' @examples
#' ZCpermit_report(2013)
#' 
ZCpermit_report=function(year,dir=NULL)
{
	Brand=getCalcurData("Zc","ZcBrand",dir=dir)
	Tags=getCalcurData("Zc","TagInitial",dir=dir)
	TagRecaps=getCalcurData("Zc","TagRecapture",dir=dir)
	BrandRecaps=getCalcurData("Zc","Recaptures",dir=dir)
	Samples=getCalcurData("Zc","SampleCollection",dir=dir)
# Extract only those records that were sampled 
	Brand=Brand[Brand$sample=="Y"&!is.na(Brand$sample),]
	Tags=Tags[!is.na(Tags$sample) & Tags$sample=="Y",]
	TagRecaps=TagRecaps[!is.na(TagRecaps$sampled)&TagRecaps$sampled=="Y",]
	BrandRecaps=BrandRecaps[!is.na(BrandRecaps$sampled)&BrandRecaps$sampled=="Y",]
# create PYear field which is the Permit Year. It is the year of capture date except when 
# month is between Jan and June (<7)
	Brand$PYear= as.POSIXlt(Brand$branddate)$year+1900
    Brand$Year=Brand$PYear
	Brand$Month=as.POSIXlt(Brand$branddate)$mon+1
	Brand$PYear[as.POSIXlt(Brand$branddate)$mon+1<7]=Brand$PYear[as.POSIXlt(Brand$branddate)$mon+1<7]-1
	Tags$PYear= as.POSIXlt(Tags$capturedate)$year+1900
	Tags$Year=Tags$PYear
	Tags$Month=as.POSIXlt(Tags$capturedate)$mon+1
	Tags$PYear[as.POSIXlt(Tags$capturedate)$mon+1<7]=Tags$PYear[as.POSIXlt(Tags$capturedate)$mon+1<7]-1
	TagRecaps$PYear= as.POSIXlt(TagRecaps$capturedate)$year+1900
	TagRecaps$Year=TagRecaps$PYear
	TagRecaps$Month=as.POSIXlt(TagRecaps$capturedate)$mon+1
	TagRecaps$PYear[as.POSIXlt(TagRecaps$capturedate)$mon+1<7]=TagRecaps$PYear[as.POSIXlt(TagRecaps$capturedate)$mon+1<7]-1
	BrandRecaps$PYear= as.POSIXlt(BrandRecaps$capturedate)$year+1900
	BrandRecaps$Year=BrandRecaps$PYear
	BrandRecaps$Month=as.POSIXlt(BrandRecaps$capturedate)$mon+1
	BrandRecaps$PYear[as.POSIXlt(BrandRecaps$capturedate)$mon+1<7]=BrandRecaps$PYear[as.POSIXlt(BrandRecaps$capturedate)$mon+1<7]-1
	
# Extract set of common fields from each data source
	Brand=subset(Brand,select=c("AnimalID","region","sex","Permit","PYear","Year","Month"))
	Brand$ageclass="PUP"
	Tags=subset(Tags,select=c("AnimalID","region","sex","Permit","PYear","ageclass","Year","Month"))
	Tags$ageclass=as.character(Tags$ageclass)
	TagRecaps=subset(TagRecaps,select=c("AnimalID","region","Sex","Permit","PYear","ageclass","Year","Month"))
	TagRecaps$ageclass=as.character(TagRecaps$ageclass)
	names(TagRecaps)[3]="sex"
	BrandRecaps=subset(BrandRecaps,select=c("AnimalID","region","sex","Permit","PYear","ageclass","Year","Month"))
	BrandRecaps$ageclass=as.character(BrandRecaps$ageclass)
# append all common data from each source
	SampledZc=rbind(Brand,Tags,TagRecaps,BrandRecaps)
	SampledZc=droplevels(SampledZc)
# patch some values of ageclass that are slightly different
	SampledZc$ageclass[SampledZc$ageclass=="P"]="PUP"
	SampledZc$ageclass[SampledZc$ageclass=="Pup"]="PUP"
	SampledZc$ageclass[SampledZc$ageclass=="J"]="JUV"
	SampledZc$ageclass[SampledZc$ageclass=="SAM"]="ADM"
	SampledZc$ageclass[SampledZc$ageclass=="Y"]="JUV"
	SampledZc$ageclass[SampledZc$ageclass=="SMADM"]="ADM"
# Merge SampledZc with Samples (SampleCollection table)
	ZcSamples=merge(SampledZc,Samples,by="AnimalID")
# Create a giant array by region,sex,ageclass,AnimalId and sampletype for a specific Permit Year
	py_sample=droplevels(ZcSamples[ZcSamples$PYear==year,])
	if(nrow(py_sample)==0) stop("\nNo samples found\n")
	# it now loops over unique set of years and months and puts into a list to return
    ym=(unique(subset(py_sample,select=c("Year","Month"))))
    ym=ym[order(ym$Year,ym$Month),,drop=FALSE]
	permit_list=NULL
    for(p in 1:nrow(ym))	
	{
		sample_table=with(py_sample[py_sample$Year==ym$Year[p]&py_sample$Month==ym$Month[p],,drop=FALSE],
				             table(region,sex,ageclass,AnimalID,sampletype))
		sample_table[sample_table>1]=1
#        some nasty code to create the needed table
		charray=array(NA, dim(sample_table)[1:4])
		for(i in 1:dim(sample_table)[1])
			for(j in 1:dim(sample_table)[2])
				for(k in 1:dim(sample_table)[3])
					for(l in 1:dim(sample_table)[4])
						charray[i,j,k,l]=paste(dimnames(sample_table)[5]$sampletype[sample_table[i,j,k,l,]>0],collapse=",") 
		permit_table=NULL
		for(i in 1:dim(charray)[1])
		{
			for(j in 1:dim(charray)[2])
				for(k in 1:dim(charray)[3])
				{
					xx=charray[i,j,k,][charray[i,j,k,]!=""]
					if(length(xx)>0)
					{
						xz=table(charray[i,j,k,][charray[i,j,k,]!=""])
						nl=length(xz)
						permit_table=rbind(permit_table,data.frame(region=rep(dimnames(sample_table)[1]$region[i],nl),sex=rep(dimnames(sample_table)[2]$sex[j],nl),ageclass=rep(dimnames(sample_table)[3]$ageclass[k],nl),samples=names(xz),frequency=as.vector(xz),stringsAsFactors=FALSE))
					}
				}	
		}
		permit_list[[p]]=permit_table
	}
	names(permit_list)=paste(ym$Year,ym$Month,sep="-")
	permit_list
}	


