# Plot results from bcftools ROH analysis
# It takes takes one table as an argument and plots several aspects of ROH
# just grep regions from the table: grep RG mocker_to_toxo_snps_fin.dphetmask.maxdp250.mindp8.vcf.gz.roh | head > ROH.rg

#install.packages("ggrepel", lib="/storage/brno3-cerit/home/vlkofly/Rpackages/", repos="https://cran.wu.ac.at/")
#library(devtools,lib.loc="/storage/brno3-cerit/home/vlkofly/Rpackages/")
#install.packages("RColorBrewer", lib="/storage/brno3-cerit/home/vlkofly/Rpackages/", repos="https://cran.wu.ac.at/")
#install_github("slowkow/ggrepel", args = c('--library="/storage/brno3-cerit/home/vlkofly/Rpackages/"'))
#install.packages("Rcpp", lib = "/storage/brno3-cerit/home/vlkofly/Rpackages/",repos="https://cran.wu.ac.at/",dependencies=T)
###
tot_length_ref <- 887012650 ## total length of the reference on which the snpset is based again have to get callable sites - depthmask
args <- commandArgs(TRUE)
inp <- args[1]
print(inp)

library(ggplot2,lib.loc="/storage/brno3-cerit/home/vlkofly/Rpackages/")
library(reshape2,lib.loc="/storage/brno3-cerit/home/vlkofly/Rpackages/")
library(ggrepel,lib.loc="/storage/brno3-cerit/home/vlkofly/Rpackages/")
library(RColorBrewer,lib.loc="/storage/brno3-cerit/home/vlkofly/Rpackages/")

roh.reg <- read.delim(inp, header=T) # supply the source table
head (roh.reg)
summary(roh.reg)
# filter good quality ROH with phred score > 30, only 10 regions
roh.reg<-roh.reg[roh.reg$X.8.Quality..average.fwd.bwd.phred.score. > 30,]
summary(roh.reg)
####per sample summary
samples<-(levels(roh.reg$X.2.Sample))
aver.length.roh<-aggregate(roh.reg$X.6.Length..bp., list(roh.reg$X.2.Sample), function(x) c(mn = mean(x), sd = sd(x), max = max(x), sum = sum(x) ))# average length of ROH per sample
averchr.length.roh<-aggregate(roh.reg$X.6.Length..bp., list(roh.reg$X.3.Chromosome), function(x) c(mn = mean(x), sd = sd(x), max = max(x)))# average length of ROH per chrom
summary(averchr.length.roh)
write.table(aver.length.roh,paste("rohlen_summary_persamples",inp,".tsv",sep=""),row.names=F,quote=F,sep="\t")
###plot distribution of ROH per samples box & whiskers persample
# do locally


#function for plotting NROH SROH based on 10.1038/nrg.2017.109 Ceballos 2018
nroh_sroh<-function(minlen,df) {  # SROH/individual >0.5M, this is sum of length of ROH longer then 0.5M
  dfrest<-df[df$X.6.Length..bp. > minlen,]
  sroh<-aggregate(dfrest$X.6.Length..bp., list(dfrest$X.2.Sample), sum)
  nroh<-aggregate(dfrest$X.6.Length..bp., list(dfrest$X.2.Sample), length)
  Froh<-sroh$x/tot_length_ref
  snroh<-cbind("sroh" = sroh, "nroh" = nroh$x, "Froh" = Froh, "length.category" = rep(minlen,length(Froh)))
  return(snroh)
  }

nroh_sroh_range<-function(minlen,maxlen,df) { # function to plot distribution of ROH in ranges lets say from 0.5M - 1M
dfrest<-df[df$X.6.Length..bp. > minlen & df$X.6.Length..bp. < maxlen,]
sroh<-aggregate(dfrest$X.6.Length..bp., list(dfrest$X.2.Sample), sum)
nroh<-aggregate(dfrest$X.6.Length..bp., list(dfrest$X.2.Sample), length)
Froh<-sroh$x/tot_length_ref
snroh<-cbind("sroh" = sroh, "nroh" = nroh$x, "Froh" = Froh, "length.category" = rep(paste(minlen,maxlen,sep="-"),length(Froh)))
	      return(snroh)

}


sord<-c("C_254","C_256","C_270","D_001","D_002","D_013","W_001","W_009","W_025","F_272","F_276","F_280","P_001","P_016","P_017","M_182","M_192","M_208","S_001","S_005","S_140","I_364","I_370","I_380","L_001","L_002","L_003")
rohlen<-c(1e5,5e5,1e6)
rohanal<-lapply(rohlen,function(x) nroh_sroh(x,roh.reg) )
rohanal<-do.call("rbind",rohanal)
rohanal$sroh.Group.1 <- factor (rohanal$sroh.Group.1, levels=sord) # reorder levels

write.table(rohanal,"nroh_distribution_longer.tsv",row.names=F,quote=F)

levels(rohanal$sroh.Group.1)
th<-theme(axis.text.x=element_text(size=15),
          axis.title = element_text(size=23),
          axis.text.y=element_text(size=15),plot.title=element_text(size=25))

col<-brewer.pal(9,"Set1")
col<-rep(col,each=3)

svg("nroh_distribution_longer.svg")
ggplot(rohanal,aes(x=length.category,y=nroh,fill=sroh.Group.1))+
  geom_bar(stat="identity",position = "dodge")+
  scale_x_discrete(limits=rohlen,labels=c(">0.1Mbp",">0.5Mbp",">1Mbp"))+
  scale_fill_manual(values=col) +
  xlab("ROH length category")+
  ylab("Number of ROH within category")+
  th
dev.off()

rohlen<-c(1e6,1.5e6,2e6)
rohanal<-lapply(rohlen,function(x) nroh_sroh(x,roh.reg))
rohanal<-do.call("rbind",rohanal)
rohanal$sroh.Group.1 <- factor (rohanal$sroh.Group.1, levels=sord) # reorder levels
write.table(rohanal,"nroh_distribution_shorter.tsv",row.names=F,quote=F)


svg("nroh_distribution_shorter.svg")
ggplot(rohanal,aes(x=length.category,y=nroh,fill=sroh.Group.1))+
  geom_bar(stat="identity",position = "dodge")+
  scale_fill_manual(values=col)+
  scale_x_discrete(limits=rohlen,labels=c(">1Mbp",">1.5Mbp",">2Mbp"))+
  xlab("ROH length category")+
  ylab("Number of ROH within category")+
  th
dev.off()

#so now print the range stats
rohanal<-rbind(nroh_sroh_range(1e5,5e5,roh.reg),nroh_sroh_range(5e5,1e6,roh.reg), nroh_sroh_range(1e6,2e6,roh.reg),nroh_sroh_range(2e6,800e6,roh.reg))
rohanal$sroh.Group.1 <- factor (rohanal$sroh.Group.1, levels=sord) # reorder levels
write.table(rohanal,"nroh_distribution_range.tsv",row.names=F,quote=F)

svg("nroh_distribution_range_wider.svg")
ggplot(rohanal,aes(x=length.category,y=nroh,fill=sroh.Group.1))+
geom_bar(stat="identity",position = "dodge")+
scale_fill_manual(values=col)+
scale_x_discrete(labels=c("0.1-0.5","0.5-1","1-2",">2"))+
xlab("ROH length category in Mbp")+
ylab("Number of ROH within category")+
th
dev.off()



roh0.5m<-nroh_sroh(5e5,roh.reg)

names(roh0.5m)
svg("ROHbcftools.nrohsroh0.5.svg")
ggplot(roh0.5m,aes(sroh.x/1e6,nroh,label = sroh.Group.1)) +
  geom_point()+
  geom_text_repel(size=5)+
  scale_x_continuous(limits=c(0,800))+
  #geom_label_repel(data=roh0.5m[roh0.5m$Froh > 0.18,],aes(x=580,label=round(Froh,2)),size=3)+
  xlab("Total ROH length in Mb (sroh) ")+
  ylab("ROH count (nroh)")+
  ggtitle("ROH distribution >0.5Mb category")+
  th

dev.off()
svg("ROH_zoom.nrohsroh0.5.svg")
ggplot(roh0.5m[roh0.5m$nroh < 100,],aes(sroh.x/1e6,nroh,label = sroh.Group.1)) +
  geom_point()+
  geom_text_repel(size=5)+
  xlab("Total ROH length in Mb (sroh) ")+
  ylab("ROH count (nroh)")+
  ggtitle("ROH distribution >0.5Mb category")+
  th
#geom_label_repel(data=roh0.5m[roh0.5m$Froh > 0.18,],aes(x=580,label=round(Froh,2)),size=3)+

dev.off()



# plot ROH some chromosome

plotroh<- function(chrom){
  roh.scaff437<-roh.reg[roh.reg$X.3.Chromosome == chrom,]
  samord<-sort(as.character(unique(roh.scaff437$X.2.Sample))) # play with this
  print(samord)
  samord<-samord[c(1:6,22:24,7:9,13:18,19:21,10:12)] # order samples have to change for every batch
 # snum<-as.numeric(gsub("-.*","",as.character(roh.scaff437$X.2.Sample))) # change sample name into integer
 # roh.scaff437<-cbind(roh.scaff437,snum)
  g<-ggplot()+
    geom_segment(data=roh.scaff437,mapping=aes(x=X.4.Start,xend=X.5.End,y=X.2.Sample,yend=X.2.Sample),size=2)+
    scale_y_discrete(limits=samord)+
    ylab('')+
    xlab("")+
    ggtitle (paste("Runs of homozygosity:",chrom))+
    theme(panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5,size=20),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15))
  print(g)
}

save.image("ROH.Rdata")
#plotroh("scaffold437")
chroms<-c("Mmel_Chr1","Mmel_Chr2","Mmel_Chr3","Mmel_Chr4","Mmel_Chr10","Mmel_Chr22","Mmel_Chr28","Mmel_ChrZ")
#chroms<-c("scaffold437","scaffold46")

pdf("LROHbcftools.pdf",width = 40)
lapply(chroms, function(x) plotroh(x))
dev.off()

