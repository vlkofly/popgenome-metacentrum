# load data specify table in the argument
# run under version of R-3.4.3-gcc

####threads###
nt =1 

# install SNPrelate
#source("https://bioconductor.org/biocLite.R")
#biocLite("SNPRelate", lib = "/storage/brno3-cerit/home/janav/Rpackages/")
#biocLite("gdsfmt", lib = "/storage/brno3-cerit/home/janav/Rpackages/")
#biocLite("GWASTools", lib = "/storage/brno3-cerit/home/vlkofly/Rpackages/")
#install.packages("devtools", lib = "/storage/brno3-cerit/home/janav/Rpackages/",repos="https://cran.wu.ac.at/")
#install.packages("Rcpp", lib = "/storage/brno3-cerit/home/janav/Rpackages/",repos="https://cran.wu.ac.at/",dependencies=T)
#install.packages("ggplot2", lib = "/storage/brno3-cerit/home/janav/Rpackages/",repos="https://cran.wu.ac.at/",dependencies=T)

#library(devtools,lib.loc="/storage/brno3-cerit/home/janav/Rpackages/")
#install_github("slowkow/ggrepel", args = c('--library="/storage/brno3-cerit/home/janav/Rpackages/"'))


library(gdsfmt,lib.loc="/storage/brno3-cerit/home/janav/Rpackages/")
library(SNPRelate,lib.loc="/storage/brno3-cerit/home/janav/Rpackages/")
library(ggplot2,lib.loc="/storage/brno3-cerit/home/janav/Rpackages/")
library(ggrepel,lib.loc="/storage/brno3-cerit/home/janav/Rpackages/")
library(GWASTools,lib.loc="/storage/brno3-cerit/home/vlkofly/Rpackages/")

args <- commandArgs(TRUE)

inputdata <- args[1]

print(inputdata)
gdsfile = gsub(".vcf",".gds",inputdata)
basename = gsub(".vcf","",inputdata)
####import data from vcf####
snpgdsVCF2GDS(inputdata,gdsfile,ignore.chr.prefix="Mmel_Chr")# I guess it needs extracted vcf modified for melanotis genome
snpgdsSummary(gdsfile)
genofile = snpgdsOpen(gdsfile)

####select unlinked loci######
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2, maf=0.05, autosome.only=FALSE, missing.rate = 0.2,  num.thread = nt) # maf filtering added August 2018 according Jacquelin work 0.1 or 0.05
snp.id=unlist(snpset)
Nsnp<-length(snp.id)
print(Nsnp)





######pca#######

pca <- snpgdsPCA(genofile, snp.id=snp.id, num.thread= nt, autosome.only=FALSE)

pc.percent <- pca$varprop*100
pc.percent <-(round(pc.percent, 2))
pc1name <-paste("PC1",pc.percent[1])
pc2name <-paste("PC2",pc.percent[2])

# make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)

svg(paste(basename,"PCA.svg",sep="")) #export to svg
ggplot(data=tab, aes(x = EV1, y = EV2, label = sample.id))+ #colour = pop.tlr
  geom_point()+
  labs(x=pc1name,y=pc2name) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_text_repel(colour = "black", alpha = 0.9, size = 3 ) +
  ggtitle(paste("pca based on",Nsnp, "unlinked snps")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


#####calculate identity by state between samples######
dissMatrix  =  snpgdsIBS(genofile , sample.id=NULL,snp.id=snp.id, autosome.only=FALSE, 
    remove.monosnp=TRUE,  maf=NaN, missing.rate=NaN, num.thread=nt, verbose=TRUE)

snpHCluster =  snpgdsHCluster(dissMatrix, sample.id=NULL, need.mat=TRUE, hang=0.01)

cutTree = snpgdsCutTree(snpHCluster, z.threshold=15, outlier.n=5, n.perm = 5000, samp.group=NULL, 
    col.outlier="red", col.list=NULL, pch.outlier=4, pch.list=NULL,label.H=FALSE, label.Z=TRUE, 
    verbose=TRUE)

svg(paste(basename,"ibstree.svg",sep=""))
snpgdsDrawTree(cutTree, main = paste("Relatedness IBS based on",Nsnp, "unlinked snps"),edgePar=list(col=rgb(0.5,0.5,0.5,0.75),t.col="black"),
    y.label.kinship=T,leaflab="perpendicular")
dev.off()

#####calculate disst between samples######
dissMatrix.d  =  snpgdsDiss(genofile , sample.id=NULL,snp.id=snp.id, autosome.only=FALSE,
    remove.monosnp=TRUE,  maf=NaN, missing.rate=NaN, num.thread=nt, verbose=TRUE)

snpHCluster.d =  snpgdsHCluster(dissMatrix.d, sample.id=NULL, need.mat=TRUE, hang=0.01)

cutTree = snpgdsCutTree(snpHCluster.d, z.threshold=15, outlier.n=5, n.perm = 5000, samp.group=NULL,
    col.outlier="red", col.list=NULL, pch.outlier=4, pch.list=NULL,label.H=FALSE, label.Z=TRUE,
    verbose=TRUE)

svg(paste(basename,"disstree.svg",sep=""))
snpgdsDrawTree(cutTree, main = paste("Relatedness diss based on",Nsnp, "unlinked snps"),edgePar=list(col=rgb(0.5,0.5,0.5,0.75),t.col="black"),
    y.label.kinship=T,leaflab="perpendicular")
dev.off()

####individual inbreeding######
sample.id = read.gdsn(index.gdsn(genofile,"sample.id"))# new


indinb = snpgdsIndInb(genofile, snp.id=snp.id, autosome.only=FALSE,
    remove.monosnp=TRUE,  maf=NaN, missing.rate=NaN, verbose=TRUE, method = "mom.visscher")

print(sample.id)
print(indinb$inbreeding)
indinb = cbind(sample.id,indinb$inbreeding)
# determine how to output it
write.table(indinb, (paste(basename,"indinb.tsv",sep="")), quote=F, sep="\t", row.names=F)


####identity by descent####

ibd  = snpgdsIBDMoM(genofile, sample.id=NULL, snp.id=snp.id, autosome.only=FALSE,
    remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, allele.freq=NULL,
    kinship=TRUE, kinship.constraint=FALSE, num.thread=nt, verbose=TRUE)

ibd.coeff = snpgdsIBDSelection(ibd)

svg(paste(basename,"relatednessMOM.svg",sep=""))
plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
    xlab="k0", ylab="k1", main="Relatedness (MoM)")
lines(c(0,1), c(1,0), col="red", lty=2)
dev.off()

write.table(ibd.coeff,paste(basename,"relatednessIBD_MoM.tsv", sep=""), quote=F, sep="\t",row.names=F)

####ibd by king####
ibdk  = snpgdsIBDKING(genofile, sample.id=NULL, snp.id=snp.id, autosome.only=FALSE,
    remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    num.thread=nt,type="KING-homo", verbose=TRUE)

ibdk.coeff = snpgdsIBDSelection(ibdk)

write.table(ibdk.coeff,paste(basename,"relatednessIBD_KING_homo.tsv", sep=""), quote=F, sep="\t",row.names=F)

save.image(paste(basename,"_snprelate.RData",sep=""))


# subset the gds file for the snpset
# snp.id goes from one
# snp.position is the thing you need to output

#snpgdsClose(genofile)
#gdsSubset(gdsfile,"neutral.unlinked.gds",snp.include=snp.id)

snpgdsGDS2BED(genofile,"unlinked.neutral.bed",snp.id=snp.id) # wow this was really simple to output it in plink



