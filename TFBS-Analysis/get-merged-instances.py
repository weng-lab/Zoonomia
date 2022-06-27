import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logomaker
import json
import re
from multiprocessing import Pool
import random
import string
from subprocess import run, Popen, PIPE
import os
import sys

def get_random_bed():
    letters = string.ascii_lowercase
    return(''.join(random.choice(letters) for i in range(10)) + ".bed")

def get_information_content(x):
    ic = x * np.log2((x + .001) / (.25 + .001))
    if ic > 0:
        return(ic)
    else:
        return(0.0)
    
def trim_ppm(ppm, min_info=0.25):
    w = ppm.shape[0]
    info = np.zeros(w)
    for i in range(w):
        for j in range(4):
            info[i] += ppm[i,j] * np.log2((ppm[i,j] + .001) / 0.25)
    
    start_index = 0
    stop_index = w
    for i in range(w):
        if info[i] < min_info:
            start_index += 1
        else:
            break
            
    for i in range(w):
        if info[w-i-1] < 0.25:
            stop_index -= 1
        else:
            break
           
    if np.max(info) < 0.25:
        return(ppm)
    else:
        return(ppm[start_index:stop_index,:])
    
from scipy.stats import pearsonr
def score(m1, m2):
    
    if np.isnan(np.sum(m1)) or np.isnan(np.sum(m1)):
        return(-9999)
    
    m1 = np.nan_to_num(m1, nan=0.25)
    m2 = np.nan_to_num(m2, nan=0.25)
    
    w1 = m1.shape[0]
    w2 = m2.shape[0]
    
    if w1 < w2:
        query = m1
        template = m2
    else:
        query = m2
        template = m1
        
        
    query_w = query.shape[0]
    template_w = template.shape[0]
    query_rc = query[::-1,::-1]
    
    pccs = []
    padded_template = 0.25 * np.ones((2*query_w + template_w, 4))
    padded_template[query_w:(query_w + template_w)] = template
    
    max_pcc = -9999
    for i in range(1,padded_template.shape[0] - query_w - 1): 
        pcc_for = pearsonr(query.flatten(), padded_template[i:(i+query_w),:].flatten())[0]
        pcc_rc = pearsonr(query_rc.flatten(), padded_template[i:(i+query_w),:].flatten())[0]
        
        if np.max((pcc_for, pcc_rc)) > max_pcc:
            max_pcc = np.max((pcc_for, pcc_rc))
            
    return(max_pcc)

from scipy.stats import pearsonr
def score(m1, m2):
    
    if np.isnan(np.sum(m1)) or np.isnan(np.sum(m1)):
        return(-9999)
    
    m1 = np.nan_to_num(m1, nan=0.25)
    m2 = np.nan_to_num(m2, nan=0.25)
    
    w1 = m1.shape[0]
    w2 = m2.shape[0]
    
    if w1 < w2:
        query = m1
        template = m2
    else:
        query = m2
        template = m1
        
        
    query_w = query.shape[0]
    template_w = template.shape[0]
    query_rc = query[::-1,::-1]
    
    pccs = []
    padded_template = 0.25 * np.ones((2*query_w + template_w, 4))
    padded_template[query_w:(query_w + template_w)] = template
    
    max_pcc = -9999
    for i in range(1,padded_template.shape[0] - query_w - 1): 
        pcc_for = pearsonr(query.flatten(), padded_template[i:(i+query_w),:].flatten())[0]
        pcc_rc = pearsonr(query_rc.flatten(), padded_template[i:(i+query_w),:].flatten())[0]
        
        if np.max((pcc_for, pcc_rc)) > max_pcc:
            max_pcc = np.max((pcc_for, pcc_rc))
            
    return(max_pcc)

metadata = {}
with open("/data/zusers/andrewsg/GTRD/GTRD-metadata.txt") as f:
    for line in f:
        split = line.strip().split("\t")
        acc, tf, cell = split
        metadata[acc] = {"tf" : tf, "cell" : cell}
GTRD_tfs = list(set([metadata[x]["tf"] for x in metadata]))

remap_metadata = {}
with open("/data/zusers/andrewsg/ReMap/ReMap-metadata.txt") as f:
    for line in f:
        split = line.strip().split("\t")
        remap_metadata[split[0]] = split[1]
        
def get_instances(x):
    tf, acc, seed = x
    try:
        if "ZNF" in tf or "ZBTB" in tf:
            with open("./ZNF/1-ZMotif/" + acc + "/" + acc + ".json") as f:
                data = json.load(f)
        else:
            with open("./1-ZMotif/" + acc + "/" + acc + ".json") as f:
                data = json.load(f)
    except:
        return(None)
    
    for i, x in enumerate(data):
        ppm = data[x]["ppm"]
        ppm = np.array(ppm)
        ppm = trim_ppm(ppm)
        n = data[x]["n_sites"]
        pcc = score(ppm, seed)
        if pcc > 0.9:
            print(pcc)
            if "ZNF" in tf or "ZBTB" in tf:
                bed = "./ZNF/1-ZMotif/" + acc + "/" + acc + ".bed"
            else:
                bed = "./1-ZMotif/" + acc + "/" + acc + ".bed"
              
            
            tmp_bed = get_random_bed()
            with open(bed) as f, open("./2-Instances/" + tf + "/" + tmp_bed, "w") as g:
                for line in f:
                    split = line.strip().split("\t")
                    motif_id = split[3].split("_")[0] + "_" + split[3].split("_")[1]
                    if x == motif_id:
                        chrom, start, stop = split[:3]
                        print(chrom, start, stop, sep="\t", file=g)
            return(tmp_bed)
    return(None)

tf = sys.argv[1]
print("Obtaining seed")
with open("seeds-refined.txt") as f:
    for line in f:
        split = line.strip().split("\t")
        if split[0] == tf:
            motif_id = split[1]
            break
            
acc = motif_id.split("_")[0]
if "ZNF" in tf or "ZBTB" in tf:
    with open("./ZNF/1-ZMotif/" + acc + "/" + acc + ".json") as f:
        data = json.load(f)
else:
    with open("./1-ZMotif/" + acc + "/" + acc + ".json") as f:
        data = json.load(f)
        
seed = trim_ppm(np.array(data[motif_id]["ppm"]))
if "ZNF" in tf or "ZBTB" in tf:
    acc_list = [x for x in metadata if metadata[x]["tf"] == tf] + [x for x in remap_metadata if remap_metadata[x] == tf]
else:
    acc_list = [x for x in metadata if metadata[x]["tf"] == tf]
print(acc_list)
args = [[tf, x, seed] for x in acc_list]
with Pool(24) as p:
    beds = p.map(get_instances, args)
    
work_dir = "/data/zusers/andrewsg/zoonomia/2-Instances/" + tf + "/"
beds = [work_dir + x for x in beds if x is not None]
print(beds)
cmd = ["cat"] + beds
with open(work_dir + tf + ".instances.bed", "w") as f:
    run(cmd, stdout=f)
  
for x in beds:
    os.remove(x)

sort = Popen(["sort", "-k1,1", "-k2,2n", work_dir + tf + ".instances.bed"], stdout=PIPE)
awk = Popen(['awk', 'BEGIN{FS=OFS="\t"} $2-100 > 0 {print $1,$2-100,$3+100}'],
            stdin=sort.stdout, stdout=PIPE)
with open(work_dir + tf + ".final.bed", "w") as f:
    run(["bedtools", "merge", "-i", "-"], stdin=awk.stdout, stdout=f)
    
os.remove(work_dir + tf + ".instances.bed" )