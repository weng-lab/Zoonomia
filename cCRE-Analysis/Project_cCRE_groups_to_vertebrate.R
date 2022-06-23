# -- Kaili Fan
# fankaili.bio@gmail.com
# Weng Lab, UMass Chan Med
# Jun, 2022

# This script is for making heatscatter for cCRE in vertebrate and project groups.

library(ggplot2)
library(umap)

args = commandArgs(trailingOnly=TRUE)
vertebrate_alignment_matrix = args[1]
mammal_alignment_matrix = args[1]
species_list = args[2]
mammal_order = args[3]


######################
# scatter - alignment for each cCRE
######################
dat = read.table(vertebrate_alignment_matrix, header = TRUE)

df <- data.frame(y = dat$align90, x = dat$align10,
                 d = densCols(dat$align90, dat$align10, colramp = colorRampPalette(rev(brewer.pal(11,"Spectral")))))
ggplot(df, aes(x=x, y=y, color=d)) +
  geom_point() + scale_color_identity() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border=element_rect(colour = "black", fill=NA),
        legend.position="none") +
  ylab("#species with ≥90% alignment")+xlab("#species with ≤10% alignment") +
  coord_cartesian(xlim=c(0,100),ylim=c(0,100))

ggsave("GRCh38-ccREs_phyloP100_alignment_scatter_rmPAR.png", width=5)


######################
# project groups to 241-mammal
######################
dat0 = read.table(vertebrate_alignment_matrix, header = TRUE)

df0 = dat0[dat0$align90>=61,]$id
length(df0)
length(df0)/nrow(dat0)
df0b = dat0[dat0$align90>=25 & dat0$align10<38,]$id
length(df0b)
length(df0b)/nrow(dat0)
df1 = dat0[dat0$align90>=25 & dat0$align10>=38,]$id
length(df1)
length(df1)/nrow(dat0)
df2 = dat0[dat0$align90>=5 & dat0$align90<=15 & dat0$align10<=70,]$id
length(df2)
length(df2)/nrow(dat0)
df3 = dat0[dat0$align10>=77,]$id
length(df3)
length(df3)/nrow(dat0)

dat1 = read.table(mammal_alignment_matrix, header = TRUE)
rownames(dat1) = dat1$id

#------------------------
dat = dat1[df1,]

df <- data.frame(y = dat$align90, x = dat$align10,
                 d = densCols(dat$align90, dat$align10, colramp = colorRampPalette(rev(brewer.pal(11,"Spectral")))))
ggplot(df, aes(x=x, y=y, color=d)) +
  geom_point() + scale_color_identity() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border=element_rect(colour = "black", fill=NA),
        legend.position="none") +
  ylab("#species with ≥90% alignment")+xlab("#species with ≤10% alignment") +
  coord_cartesian(xlim=c(0,240),ylim=c(0,240)) +
  labs(title=paste("G1 cCREs\n(N=",as.character(length(df1)),")",sep=""))

ggsave("GRCh38-G1_phyloP240_alignment_scatter_rmPAR.png", width=5)


#------------------------
dat = dat1[df2,]

df <- data.frame(y = dat$align90, x = dat$align10,
                 d = densCols(dat$align90, dat$align10, colramp = colorRampPalette(rev(brewer.pal(11,"Spectral")))))
ggplot(df, aes(x=x, y=y, color=d)) +
  geom_point() + scale_color_identity() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border=element_rect(colour = "black", fill=NA),
        legend.position="none") +
  ylab("#species with ≥90% alignment")+xlab("#species with ≤10% alignment") +
  coord_cartesian(xlim=c(0,240),ylim=c(0,240)) +
  labs(title=paste("G2' cCREs\n(N=",as.character(length(df2)),")",sep=""))

ggsave("GRCh38-G2'_phyloP240_alignment_scatter_rmPAR.png", width=5)

#------------------------
dat = dat1[df3,]

df <- data.frame(y = dat$align90, x = dat$align10,
                 d = densCols(dat$align90, dat$align10, colramp = colorRampPalette(rev(brewer.pal(11,"Spectral")))))
ggplot(df, aes(x=x, y=y, color=d)) +
  geom_point() + scale_color_identity() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border=element_rect(colour = "black", fill=NA),
        legend.position="none") +
  ylab("#species with ≥90% alignment")+xlab("#species with ≤10% alignment") +
  coord_cartesian(xlim=c(0,240),ylim=c(0,240)) +
  labs(title=paste("G3 cCREs\n(N=",as.character(length(df3)),")",sep=""))

ggsave("GRCh38-G3_phyloP240_alignment_scatter_rmPAR.png", width=5)

#------------------------
dat = dat1[df0,]

df <- data.frame(y = dat$align90, x = dat$align10,
                 d = densCols(dat$align90, dat$align10, colramp = colorRampPalette(rev(brewer.pal(11,"Spectral")))))
ggplot(df, aes(x=x, y=y, color=d)) +
  geom_point() + scale_color_identity() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border=element_rect(colour = "black", fill=NA),
        legend.position="none") +
  ylab("#species with ≥90% alignment")+xlab("#species with ≤10% alignment") +
  coord_cartesian(xlim=c(0,240),ylim=c(0,240)) +
  labs(title=paste("G0 cCREs\n(N=",as.character(length(df4)),")",sep=""))

ggsave("GRCh38-G0_phyloP240_alignment_scatter_rmPAR.png", width=5)

#------------------------

dat = dat1[df0b,]

df <- data.frame(y = dat$align90, x = dat$align10,
                 d = densCols(dat$align90, dat$align10, colramp = colorRampPalette(rev(brewer.pal(11,"Spectral")))))
ggplot(df, aes(x=x, y=y, color=d)) +
  geom_point() + scale_color_identity() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border=element_rect(colour = "black", fill=NA),
        legend.position="none") +
  ylab("#species with ≥90% alignment")+xlab("#species with ≤10% alignment") +
  coord_cartesian(xlim=c(0,240),ylim=c(0,240)) +
  labs(title=paste("G0' cCREs\n(N=",as.character(length(df4b)),")",sep=""))

ggsave("GRCh38-G0'_phyloP240_alignment_scatter_rmPAR.png", width=5)
