
# -- Kaili Fan
# fankaili.bio@gmail.com
# Weng Lab, UMass Chan Med
# Jun, 2022

# This script is for making aggregation plot for given rPeak signal matrix.

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
narrow_rPeak_matrix = args[1]
broad_rPeak_matrix = args[2]

####################
dat1 = read.table(narrow_rPeak_matrix, row.names = 1)
dat2 = read.table(broad_rPeak_matrix, row.names = 1)

d1 = data.frame(apply(dat1,2,mean))
d1$peak = "narrow"
d1$loci = 1:500
#
d2 = data.frame(apply(dat2,2,mean))
d2$peak = "broad"
d2$loci = 1:500

colnames(d1) = colnames(d2) = c("score", "peak","loci")
d = data.frame(rbind(d1,d2))
d$peak = factor(d$peak, levels=c("narrow","broad"))

ggplot(d, aes(x=loci, y=score,col=peak)) +
  geom_smooth(method="loess",span=0.1, se=FALSE) +
  theme_classic() +
  xlab("") + ylab("phyloP score") +
  scale_x_continuous(breaks=c(0,250,500),
                     labels=c("-250bp","center","250bp")) +
  geom_vline(xintercept = 250, linetype="dashed", col="gray") +
  scale_color_manual(values=c("#e6ab02","#66a61e"))

ggsave(output_fig, width=5)
