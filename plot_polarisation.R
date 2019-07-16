# this R script takes variant table from GATK command VariantsToTable and generates graphs of distribution of different variant parameters

# load data specify table in the argument

#install.packages("labeling", lib="/storage/brno3-cerit/home/vlkofly/Rpackages/", repos="https://cran.wu.ac.at/")
#install.packages("grid", lib="/storage/brno3-cerit/home/vlkofly/Rpackages/", repos="https://cran.wu.ac.at/")
#install.packages("data.table", lib="/storage/brno3-cerit/home/vlkofly/Rpackages/", repos="https://cran.wu.ac.at/")

###the script has been changed to parse non-variants
### second argument added if novariants analyse then supply "novar"

library(labeling,lib.loc="/storage/brno3-cerit/home/vlkofly/Rpackages/")
library(ggplot2,lib.loc="/storage/brno3-cerit/home/vlkofly/Rpackages/")
library(gridExtra,lib.loc="/storage/brno3-cerit/home/vlkofly/Rpackages/")
library(grid,lib.loc="/storage/brno3-cerit/home/vlkofly/Rpackages/")
library(reshape2,lib.loc="/storage/brno3-cerit/home/vlkofly/Rpackages/")
library(data.table,lib.loc="/storage/brno3-cerit/home/vlkofly/Rpackages/")

args <- commandArgs(TRUE)

inputdata <- args[1]

vt <- read.table(inputdata,sep="\t", header=F) # supply the source table
base<-unlist(strsplit(inputdata,".",fixed=T))[1]
dim(vt)
#sapply(vt,mean)

# add header
names(vt)<-c("chr","pos_bed","pos_vcf","outgrp_bp","target_fr","out1_fr","out2_fr","out3_fr","pos_est","prob")
#names(vt)<-c("chr","pos_bed","pos_vcf","outgrp_bp","target_fr","out1_fr","out2_fr","pos_est","prob")

ntot<-dim(vt)[1]

head (vt)

suma<-summary(vt[,c("prob","outgrp_bp")])


vartoplot <-c("prob")
# how to plot the density function with different axis limit for each variable.

plotdensity<-function(var) {
	title<-"Prob major=ancestral" # develop more this script, add some stats and so on
	p<-ggplot(vt,aes_string(x=var))+geom_histogram() + theme(axis.text = element_text(size=20))+
	ggtitle(title)
}




svg(paste(base,"prob_hist.svg",sep="."))
lapply(vartoplot,plotdensity)
dev.off()


###make polarisation key
key<-vt[vt$prob < 0.4,c(1,2,3)] # the key is in bed format
write.table(key,paste(base,"polarisation.key",sep="."), quote=F,row.names=F,col.names=F,sep="\t")
cat ("length of polarisation key = ",dim(key)[1],"\n")
####print uSFS
sfs<-read.table(paste(base,".sfs.output.txt",sep=""),sep="\t",header=F)
sfs<-cbind(sfs,1:dim(sfs)[1])
names(sfs)<-c("frequency","derived_allele_count")
head (sfs)
png(paste(base,".sfs.png",sep=""), width=680, height=360 )
ggplot(sfs,aes(x=derived_allele_count,y=frequency))+geom_bar(stat="identity") +
       ggtitle(base)+
              scale_x_continuous()
dev.off()



cat("------------------------------------------------------------\n")
print (suma)
cat("------------------------------------------------------------\n")
#cat("stats: Total length of alignment",totlen, "bp in", ntot, "blocks \n")
#cat("stats: Total length of liftover", lenconserv, " bp in", numconserv, "blocks \n")
