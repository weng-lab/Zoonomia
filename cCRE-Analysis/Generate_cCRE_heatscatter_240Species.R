
# -- Kaili Fan
# fankaili.bio@gmail.com
# Weng Lab, UMass Chan Med
# Jun, 2022

# This script is for making heatscatter for cCREs.

library(ggplot2)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)
alignment_count = args[1]
alignment_count_random = args[2]

######################
# heatscatter - alignment for each cCRE
######################
dat = read.table(alignment_count, header = TRUE)

df <- data.frame(x = dat$align90, y = dat$align10,
                 d = densCols(dat$align90, dat$align10, colramp = colorRampPalette(rev(brewer.pal(11,"Spectral")))))
ggplot(df, aes(x=x, y=y, color=d)) +
  geom_point() + scale_color_identity() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border=element_rect(colour = "black", fill=NA),
        legend.position="none") +
  xlab("#species with ≥90% alignment")+ylab("#species with ≤10% alignment") +
  coord_cartesian(xlim=c(0,240),ylim=c(0,240))

ggsave("GRCh38-ccREs_alignment_scatter_rmPAR.png", width=5)

#--------------------
ccre_list = c("PLS","pELS","dELS","DNase-H3K4me3","CTCF-only")

for(i in 1:5){
  ccre = ccre_list[i]

  dat0 = dat[dat$ccre==ccre,]
  title = paste(ccre,"\n(N=",nrow(dat0),")",sep="")
  df <- data.frame(x = dat0$align90, y = dat0$align10,
                   d = densCols(dat0$align90, dat0$align10, colramp = colorRampPalette(rev(brewer.pal(11,"Spectral")))))
  ggplot(df, aes(x=x, y=y, color=d)) +
    geom_point() + scale_color_identity() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border=element_rect(colour = "black", fill=NA),
          legend.position="none") +
    xlab("#species with ≥90% alignment")+ylab("#species with ≤10% alignment") +
    coord_cartesian(xlim=c(0,240),ylim=c(0,240))

  ggsave(paste("GRCh38-",ccre,"_alignment_scatter_rmPAR.png", sep=""),width=5)
}


######################
# heatscatter - randon regions
######################

dat = read.table(alignment_count_random, header = TRUE)

df <- data.frame(x = dat$align90, y = dat$align10,
                 d = densCols(dat$align90, dat$align10, colramp = colorRampPalette(rev(brewer.pal(11,"Spectral")))))
ggplot(df, aes(x=x, y=y, color=d)) +
  geom_point() + scale_color_identity() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border=element_rect(colour = "black", fill=NA),
        legend.position="none") +
  xlab("#species with ≥90% alignment")+ylab("#species with ≤10% alignment") +
  coord_cartesian(xlim=c(0,240),ylim=c(0,240))

ggsave("GRCh38_randomRegion_alignment_scatter.png", width=5)
