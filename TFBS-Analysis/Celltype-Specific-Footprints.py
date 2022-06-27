import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyBigWig
import tqdm
from multiprocessing import Pool
import os

def get_bigwig_vals_over_bed(x, w=150):
    bed, bigwig = x
    vals = []
    with open(bed) as f:
        lines = [line.strip().split() for line in f.readlines()]
    
    n = len(lines)
    with pyBigWig.open(bigwig) as bw:
        for i in tqdm.trange(n, disable=False):
            split = lines[i]
            chrom, start, stop = split[:3]
            start, stop = int(start), int(stop)
            center = (stop + start) // 2
            try:
                vals.append(np.nan_to_num(bw.values(chrom, center - w, center + w)))
            except:
                vals.append(np.zeros(2*w))
            
    return(np.array(vals))


args = []
bed_files = ["./data/" + x for x in os.listdir("./data2/") if "bed" in x]
bigwigs = ["./data/HepG2.bigWig", "./data/A549.bigWig"]
for x in bigwigs:
    for y in bed_files:
        args.append([y,x])
        
with Pool(16) as p:
    vals = p.map(get_bigwig_vals_over_bed, args)
    
plt.rcParams['pdf.fonttype'] = 42 
plt.rcParams['svg.fonttype'] = 'none'

colors = {0 : "C2", 1 : "C0", 2 : "C1", 3 : "C6", 4 : "C5", 5 : "C7", 6 : "C4", 7 : "C3"}
fig, ax = plt.subplots(figsize=(8,8))
for i, bed in enumerate(bed_files):
    print(bed)
    y = np.mean(vals[i] ,axis=0)
    x = np.arange(300)
    label = bed.split("/")[2].strip(".bed")
    ax.plot(x,y,label=label, color=colors[i])
ax.legend()
ax.set_ylabel("HepG2 DNase (Henry)")
ax.set_xticks([0,150,300])
ax.set_xticklabels(["-150bp", "TFBS center", "+150bp"])
ax.set_title("Cell-type specific footprints in HepG2 (Henry)")

ratio = 1.0
xleft, xright = ax.get_xlim()
ybottom, ytop = ax.get_ylim()
ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)



# Hide the right and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')

plt.savefig("HepG2.pdf")

fig, ax = plt.subplots(figsize=(8,8))
for i, bed in enumerate(bed_files):
    print(bed)
    y = np.mean(vals[8+i], axis=0)
    x = np.arange(300)
    label = bed.split("/")[2].strip(".bed")
    ax.plot(x,y,label=label, color=colors[i])
ax.legend()
ax.set_ylabel("A549 DNase (Henry)")
ax.set_xticks([0,150,300])
ax.set_xticklabels(["-150bp", "TFBS center", "+150bp"])
ax.set_title("Cell-type specific footprints in A549 (Henry)")

ratio = 1.0
xleft, xright = ax.get_xlim()
ybottom, ytop = ax.get_ylim()
ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)


# Hide the right and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')

plt.savefig("A549.pdf")