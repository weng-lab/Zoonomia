import numpy as np
import pandas as pd
import os
from subprocess import PIPE, Popen, run
from multiprocessing import Pool
import sys

species_list = [x.strip(".bw") for x in os.listdir("/home/andrewsg/data/zoonomia/") if ".bw" in x]
species_list.remove("Homo_sapiens")

tf = sys.argv[1]
data = pd.read_csv("/data/zusers/andrewsg/zoonomia/7-Constrained-TFBSs-Core-Only/" + tf + ".bed",
                   usecols=[_ for _ in range(6)], sep="\t", header=None, names = ["chrom", "start", "stop", "name", "score", "strand"])

data[["chrom", "start", "stop"]].to_csv("tmp.bed", sep="\t", index=False, header=None)

def get_align_perc(species):
    print(species)
    coverage = Popen(["bedtools", "coverage",
                      "-a", "tmp.bed",
        "-b", "/home/andrewsg/data/zoonomia/BED/" + species + ".bg",
         "-sorted"], stdout=PIPE)

    p = []
    for line in coverage.stdout:
        split = str(line, encoding="utf-8").strip().split("\t")
        p.append(float(split[-1]))
    return(np.array(p))

with Pool(16) as p:
    X = p.map(get_align_perc, species_list)
    
X  = np.transpose(X)

data = pd.DataFrame(X, columns=species_list)
N1 = np.sum(1 * data.ge(.9).values, axis=1)
N2 = np.sum(1 * data.le(.1).values, axis=1)

data = pd.read_csv("/data/zusers/andrewsg/zoonomia/7-Constrained-TFBSs-Core-Only/" + tf + ".bed",
                   usecols=[_ for _ in range(6)], sep="\t", header=None, names = ["chrom", "start", "stop", "name", "score", "strand"])


data["N1"] = N1
data["N2"] = N2
data["group"] = "other"
data.loc[(data["N1"] >= 120) & (data["N2"] <= 25), "group"] = "group-1"
data.loc[(data["N1"] <= 50) & (data["N1"] >=20) & (data["N2"]  <= 120), "group"] = "group-2"
data.loc[(data["N1"] <= 50) & (data["N2"]  >= 180), "group"] = "group-3"

data[["N1", "N2", "group"]].to_csv(tf + ".N1-N2-group.txt", sep="\t", index=False, header=False)
