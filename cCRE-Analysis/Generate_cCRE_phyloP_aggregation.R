
# -- Kaili Fan
# fankaili.bio@gmail.com
# Weng Lab, UMass Chan Med
# Jun, 2022

# This script is for making aggregation plot for given cCRE phyloP matrix.

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
phyloP100_matrix = args[1]
phyloP241_matrix = args[2]
cCRE_file = args[3]
output_fig = args[4]

####################
dat1 = read.table(phyloP100_matrix, row.names = 1)
dat2 = read.table(phyloP241_matrix, row.names = 1)
df = read.table(cCRE_file)

dat1_pls = dat1[df[df$V5=="PLS",]$V4,]
d1_pls = data.frame(apply(dat1_pls,2,mean))
d1_pls$ccre = "PLS"
d1_pls$loci = 1:500
d1_pls$phyloP = "phyloP100"
#
dat1_pels = dat1[df[df$V5=="pELS",]$V4,]
d1_pels = data.frame(apply(dat1_pels,2,mean))
d1_pels$ccre = "pELS"
d1_pels$loci = 1:500
d1_pels$phyloP = "phyloP100"
#
dat1_dels = dat1[df[df$V5=="dELS",]$V4,]
d1_dels = data.frame(apply(dat1_dels,2,mean))
d1_dels$ccre = "dELS"
d1_dels$loci = 1:500
d1_dels$phyloP = "phyloP100"
#
dat1_dnase = dat1[df[df$V5=="DNase-H3K4me3",]$V4,]
d1_dnase = data.frame(apply(dat1_dnase,2,mean))
d1_dnase$ccre = "DNase-H3K4me3"
d1_dnase$loci = 1:500
d1_dnase$phyloP = "phyloP100"
#
dat1_ctcf = dat1[df[df$V5=="CTCF-only",]$V4,]
d1_ctcf = data.frame(apply(dat1_ctcf,2,mean))
d1_ctcf$ori = "CTCF-only"
d1_ctcf$ccre = 1:500
d1_ctcf$phyloP = "phyloP100"
#
dat2_pls = dat2[df[df$V5=="PLS",]$V4,]
d2_pls = data.frame(apply(dat2_pls,2,mean))
d2_pls$ccre = "PLS"
d2_pls$loci = 1:500
d2_pls$phyloP = "phyloP241"
#
dat2_pels = dat2[df[df$V5=="pELS",]$V4,]
d2_pels = data.frame(apply(dat2_pels,2,mean))
d2_pels$ccre = "pELS"
d2_pels$loci = 1:500
d2_pels$phyloP = "phyloP241"
#
dat2_dels = dat2[df[df$V5=="dELS",]$V4,]
d2_dels = data.frame(apply(dat2_dels,2,mean))
d2_dels$ccre = "dELS"
d2_dels$loci = 1:500
d2_dels$phyloP = "phyloP241"
#
dat2_dnase = dat2[df[df$V5=="DNase-H3K4me3",]$V4,]
d2_dnase = data.frame(apply(dat2_dnase,2,mean))
d2_dnase$ccre = "DNase-H3K4me3"
d2_dnase$loci = 1:500
d2_dnase$phyloP = "phyloP241"
#
dat2_ctcf = dat2[df[df$V5=="CTCF-only",]$V4,]
d2_ctcf = data.frame(apply(dat2_ctcf,2,mean))
d2_ctcf$ori = "CTCF-only"
d2_ctcf$ccre = 1:500
d2_ctcf$phyloP = "phyloP241"

colnames(d1_pls) = colnames(d1_pels) = colnames(d1_dels) = colnames(d1_dnase) = colnames(d1_ctcf) = c("score", "ccre","loci","phyloP")
colnames(d2_pls) = colnames(d2_pels) = colnames(d2_dels) = colnames(d2_dnase) = colnames(d2_ctcf) = c("score", "ccre","loci","phyloP")
d = data.frame(rbind(d1_pls,d1_pels,d1_dels,d1_dnase,d1_ctcf,
                     d2_pls,d2_pels,d2_dels,d2_dnase,d2_ctcf))
d$ccre = factor(d$ccre, levels=c("PLS","pELS","dELS","DNase-H3K4me3","CTCF-only"))
d$phyloP = factor(d$phyloP, levels=c("phyloP241","phyloP100"))

ggplot(d, aes(x=loci, y=score,col=ccre, linetype=phyloP)) +
  geom_smooth(method="loess",span=0.1, se=FALSE) +
  theme_classic() +
  xlab("") + ylab("phyloP score") +
  scale_x_continuous(breaks=c(0,250,500),
                     labels=c("-250bp","center","250bp")) +
  geom_vline(xintercept = 250, linetype="dashed", col="gray") +
  scale_color_manual(values=c("#FF0000","#FFA700","#FFCD00","#FFAAAA","#00b0f0")) +
  geom_hline(yintercept = 0.186608, col="#a1a1a1", linetype="solid") +
  geom_hline(yintercept = 0.0817535, col="#a1a1a1", linetype="dashed")

ggsave(output_fig, width=5)
