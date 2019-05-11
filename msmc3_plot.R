# load data specify table in the argument
#install.packages("scales", lib = "/storage/brno3-cerit/home/vlkofly/Rpackages/",repos="https://cran.wu.ac.at/")

#library(labeling,lib.loc="/storage/brno3-cerit/home/janav/Rpackages/")
library(ggplot2,lib.loc="/storage/brno3-cerit/home/vlkofly/Rpackages/")
library(scales,lib.loc="/storage/brno3-cerit/home/vlkofly/Rpackages/")


args <- commandArgs(TRUE)

inputdata <- args[1]
bs <- args[2]


print(inputdata)

mu<-4.6e-9 # mutation rate of birds taken from Smeds 2016
gen<-4.5 # Grant 2000



vt <- read.table(inputdata,sep="\t", header=T) # supply the source table
n <- sub(".txt","",inputdata)
summary(vt)

bst<-read.table(bs,sep="\t", header=F)
names(bst)<-names(vt)

vt$Ne <- (1/vt$lambda)/(2*mu) #lambda_00 change in between msmc2 and msmc
vt$LeftYears <- gen*(vt$left_time_boundary/mu)
vt$RightYears <- gen*(vt$right_time_boundary/mu)

bst$Ne <- (1/bst$lambda)/(2*mu) #lambda_00 change in between msmc2 and msmc
bst$LeftYears <- gen*(bst$left_time_boundary/mu)
bst$RightYears <- gen*(bst$right_time_boundary/mu)


summary(vt)

th<-theme(axis.text.x=element_text(size=23,hjust=1),
    axis.title = element_text(size=25),
    axis.text.y=element_text(size=23),plot.title=element_text(size=25)
    )


png(paste(n,".msmc.png",sep=""),width=800, height=480,pointsize=15)
ggplot(vt,aes(x=LeftYears,y=Ne))+
  
  geom_jitter(data=bst,aes(x=LeftYears,y=Ne),stat="identity",alpha=0.3,shape=5)+

  geom_step(stat="identity")+

  theme_bw()+
  theme(legend.title = element_blank())+
  
  xlab("Years Ago ")+
  ylab("Ne")+
  #scale_y_log10(breaks=c(1e3,1e4,1e5,1e6),labels=comma)+
  scale_y_log10(labels=comma)+
  #scale_x_log10(breaks=c(1e4,1e5,1e6,5e6),labels=comma)+
  scale_x_log10(labels=comma)+
  annotation_logticks()+
  th

dev.off()

write.table(vt,paste(n,"table_ne_years.tsv",sep=""),row.names=F,col.names=T,quote=F)
