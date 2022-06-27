from collections import Counter
import numpy as np
import pandas as pd
import random
import string
import pyBigWig
from multiprocessing import Pool
from random import randint
from time import sleep
from subprocess import run
import os

[os.remove("/tmp/andrewsg/" + x) for x in os.listdir("/tmp/andrewsg/")]

non_cancer = []
with open("/data/zusers/fankaili/ccre/dname/human_WGBS_rerun/human_WGBS_non-cancer_sample_ID.txt") as f:
    for line in f:
        non_cancer.append(line.strip().split()[0])
        
wgbs_metadata = {}
with open("/data/zusers/fankaili/ccre/dname/human_WGBS_samples.txt") as f:
    for line in f:
        split = line.strip().split("\t")
        biosample, exp = split[:2]
        wgbs_metadata[exp] = biosample
        
biosamples = [wgbs_metadata[x] for x in wgbs_metadata]
experiments = [x for x in wgbs_metadata]

cancer = [x for x in experiments if x not in non_cancer]

def get_random_bed():
    letters = string.ascii_lowercase
    tmp_dir = "/tmp/andrewsg/"
    return(tmp_dir + ''.join(random.choice(letters) for i in range(10)) + ".bed")

def get_data(x):
    sleep(randint(1,10))
    bed, exp = x
    bigwig = "/data/zusers/fankaili/ccre/dname/human_WGBS_rerun/signal/" + exp + ".bw"
    tmp_bed = get_random_bed()
    cmd = ["bigWigAverageOverBed", bigwig, bed, tmp_bed]
    run(cmd, check=True)
    vals = []
    with open(tmp_bed) as f:
        for line in f:
            split = line.strip().split("\t")
            if int(split[2]) == 0:
                vals.append(np.NaN)
            else:
                vals.append(float(split[-1]))
    os.remove(tmp_bed)
    return(np.array(vals))

tmp_bed = "/tmp/andrewsg/tmp.bed"

with open("TFBSs.not-constrained.merge.bed") as f, open(tmp_bed,  "w") as tmp:
    for i, line in enumerate(f):
        split = line.strip().split("\t")
        chrom, start, stop = split
        tfbs = "TFBS_" + str(i)
        print(chrom, start, stop, tfbs, file=tmp, sep="\t")
        
args = [[tmp_bed, x] for x in experiments]
with Pool(32) as p:
    vals = p.map(get_data, args)
    
X = np.transpose(np.array(vals))
data = pd.read_csv("TFBSs.not-constrained.merge.bed", 
                   names=["chrom", "start", "stop"],
                   sep = "\t",
                   header=None)

data["median"] = np.nanmedian(X, axis=1)
data["max"] = np.nanmax(X, axis=1)
data["min"] = np.nanmin(X, axis=1)
data["range"] = data["max"] - data["min"]

indices = [i for i,x in enumerate(experiments) if x in non_cancer]
data["non_median"] = np.nanmedian(X[:,indices], axis=1)
data["non_max"] = np.nanmax(X[:,indices], axis=1)
data["non_min"] = np.nanmin(X[:,indices], axis=1)
data["non_range"] = data["non_max"] - data["non_min"]

indices = [i for i,x in enumerate(experiments) if x in cancer]
data["cancer_median"] = np.nanmedian(X[:,indices], axis=1)
data["cancer_max"] = np.nanmax(X[:,indices], axis=1)
data["cancer_min"] = np.nanmin(X[:,indices], axis=1)
data["cancer_range"] = data["cancer_max"] - data["cancer_min"]

data[["median", "range"]].dropna().to_csv("TFBSs-Not-Constrained-All-CpG-Median-Range.txt", sep="\t", index=False, header=False)
data[["non_median", "non_range"]].dropna().to_csv("TFBSs-Not-Constrained-Non-Cancer-CpG-Median-Range.txt", sep="\t", index=False, header=False)
data[["cancer_median", "cancer_range"]].dropna().to_csv("TFBSs-Not-Constrained-Cancer-CpG-Median-Range.txt", sep="\t", index=False, header=False)
