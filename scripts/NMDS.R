library(vegan)
library(ggplot2)
cluster_freq <- read.csv("../Results/cluster_freqR.csv",header=TRUE,row.names=1)
cluster_freq <- t(cluster_freq)
cluster_freqP <- cluster_freq/rowSums(cluster_freq)
meta <- read.csv("../data/Meta.csv",header=TRUE,row.names=1)

cluster_freqP.nmds <- metaMDS(cluster_freqP)

nmds_df<-scores(cluster_freqP.nmds,display=c("sites"))

nmds_df<-cbind.data.frame(nmds_df,meta)

p<-ggplot(data=nmds_df,aes(NMDS1,NMDS2,colour=Diagnosis)) + geom_point(size=3) + theme_bw()
pdf("NMDS.pdf")
plot(p)
dev.off()

nT <- ncol(cluster_freqP)

cluster_freq1 <- cluster_freq + 1.0e-6
cluster_freq1P <- cluster_freq1/rowSums(cluster_freq1)
log_cluster_freq1P <- log(cluster_freq1P)

p <- rep(0,nT)
chi <- rep(0,nT)

meanH <- rep(0,nT)
meanCD <- rep(0,nT)

for(i in 1:nT){
  temp <- kruskal.test(log_cluster_freq1P[,i] ~ meta$Diagnosis)
  tmean <- aggregate(cluster_freqP[,i]~meta$Diagnosis, FUN=mean)
  chi[i] <- temp$statistic
  p[i] <- temp$p.value
  meanCD[i] <- tmean[[2]][1]
  meanH[i] <- tmean[[2]][2]
}


pa <- p.adjust(p, method = "BH") 
sig.df <- data.frame(chi,p,pa,meanH,meanCD)
rownames(sig.df) <- colnames(log_cluster_freq1P)

scg <- read.csv("Results/clustering_gt1000_scg.csv",header=TRUE,row.names=1)
scgR <- scg[rownames(sig.df),]
sig.df$comp <- rowSums(scgR == 1)/36

