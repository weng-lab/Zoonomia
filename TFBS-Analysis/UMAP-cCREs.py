import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyBigWig
import os
import multiprocessing
from umap import UMAP
from subprocess import PIPE, Popen, run

species_list = [x.strip(".bw") for x in os.listdir("/home/andrewsg/data/zoonomia/") if ".bw" in x]
species_list.remove("Homo_sapiens")

def get_perc_align(x):
    bed, species = x
    bigwig = "/home/andrewsg/data/zoonomia/" + species + ".bw"
    vals = []
    with open(bed) as f, pyBigWig.open(bigwig) as bw:
        for line in f:
            split = line.strip().split("\t")
            chrom, start, stop = split[:3]
            start, stop = int(start), int(stop)
            l = stop - start
            try:
                y = np.nan_to_num(bw.values(chrom, start, stop))
            except:
                vals.append(0)
                continue
            
            vals.append(np.sum(y) / l)
    return(np.array(vals))

X = np.transpose(X)

u = UMAP().fit_transform(X)

np.save("cCREs-UMAP.npy", u)