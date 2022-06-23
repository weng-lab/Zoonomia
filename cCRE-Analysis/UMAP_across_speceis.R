
# -- Kaili Fan
# fankaili.bio@gmail.com
# Weng Lab, UMass Chan Med
# Jun, 2022

# This script is for making UMAP on alignment of cCREs across species.

library(ggplot2)
library(umap)

args = commandArgs(trailingOnly=TRUE)
alignment_matrix = args[1]
species_list = args[2]
mammal_order = args[3]


######################
# UMAP
######################
dat = read.table(alignment_matrix, row.names = 1, header = 1)
dat_t = t(dat)
df = read.table(species_list,row.names = 1)
color_list=read.table(mammal_order,comment.char = "")

s=1 # seed
#--- UMAP parameter
n=30 #n_neighbors
m=0.1 #min_dist

set.seed(s)
umap_model = umap(dat_t, n_neighbors=n, min_dist=m)
d_umap = as.data.frame(umap_model$layout)
d_umap$species = rownames(d_umap)
d_umap$name = df[rownames(d_umap),]$V2
d_umap$family = df[rownames(d_umap),]$V3
d_umap$order = df[rownames(d_umap),]$V4

d_umap$order = factor(d_umap$order, levels=as.vector(color_list$V1))
umap = ggplot(d_umap, aes(x=V1, y=V2, text=paste(species,family,sep="\n"), col=order)) +
  geom_point(size=1.5) +
  theme_classic() +
  scale_color_manual(name="Order",limits=color_list$V1,
                       values=color_list$V2) +
  guides(color = guide_legend(ncol=2)) +
  xlab("UMAP-1") + ylab("UMAP-2") +
  labs(subtitle=paste("n_neighbors=",n,", min_dist=",m,", seed=",s, sep=""))
umap

ggsave("GRCh38-ccREs_alignment_species-UMAP.png", umap, width=6)
ggsave("GRCh38-ccREs_alignment_species-UMAP.pdf", umap, width=6)
