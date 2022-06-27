from matplotlib import pyplot as plt
from matplotlib import cm
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from sklearn import metrics
import random
import sys
import pyBigWig
import tqdm
import logomaker
from mpl_toolkits.axes_grid.inset_locator import inset_axes


plt.rcParams['pdf.fonttype'] = 42 
plt.rcParams['svg.fonttype'] = 'none'

def get_bigwig_vals_over_bed(x, w=100):
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
            try:
                vals.append(np.nan_to_num(bw.values(chrom, start - w, stop + w)))
            except:
                vals.append(np.zeros(2*w + (stop-start)))
            
    return(np.array(vals))

def encode_sequence(seq):
    DNA_SEQ_DICT = {
    'A' : [1, 0, 0, 0],
    'C' : [0, 1, 0, 0],
    'G' : [0, 0, 1, 0],
    'T' : [0, 0, 0, 1],
    'N' : [0.25, 0.25, 0.25, 0.25]
    }
    
    return(np.array([DNA_SEQ_DICT[nuc] for nuc in seq]))

def get_information_content(x):
    ic = x * np.log2((x + .001) / .25)
    if ic > 0:
        return(ic)
    else:
        return(0.0)
    
tf = sys.argv[1]

bed = "/data/zusers/andrewsg/zoonomia/5-TFBSs-Annotated/" + tf + ".bed"

phylop = get_bigwig_vals_over_bed((bed, "/home/andrewsg/data/zoonomia/PhyloP/241-mammalian-2020v2.bigWig"))

data = pd.read_csv(bed, 
                   names = ["chrom", "start", "stop", "name", "score", "strand", "annotation", "TSS-d", "rDHS", "te", "te_class", "div", "seq"],
                   sep="\t")

n_p = len(data[data["annotation"] == "TSS"])
n_d = len(data[data["annotation"] != "TSS"])
print(n_p, n_d)
w = len(data.loc[0,"seq"])

pfm = np.sum(data["seq"].map(encode_sequence))
ppm = pfm / np.sum(pfm, axis=1, keepdims=True)

indices = data.index[(data["annotation"] != "exonic") & (data["strand"] == "+")]
y1 = np.mean(phylop[indices], axis=0)
indices = data.index[(data["annotation"] != "exonic") & (data["strand"] == "-")]
y2 = np.mean(phylop[indices], axis=0)[::-1]
y = np.mean((y1,y2), axis=0)


indices = data.index[data["annotation"] != "exonic"]
phylop_baseline = np.mean((np.mean(phylop[indices,50:100]), np.mean(phylop[indices,100+w:100+w+50])))
conserved_x = [i for i,x in enumerate(y[100:100+w]) if x > phylop_baseline]

data["phylop-flank"] = np.mean((np.mean(phylop[:,75:100], axis=1), np.mean(phylop[:,100+w:100+w+25], axis=1)), axis=0)
data["phylop-core"] = 0

rows = data.index[data["strand"] == "+"].tolist()
cols = [100+x for x in conserved_x]
data.loc[data["strand"] == "+", "phylop-core"] = np.mean(phylop[np.ix_(rows,cols)], axis=1)

rows = data.index[data["strand"] == "-"].tolist()
cols = [100+(w-x-1) for x in conserved_x]
data.loc[data["strand"] == "-", "phylop-core"] = np.mean(phylop[np.ix_(rows,cols)], axis=1)

data["constrained-score"] = data["phylop-core"] - data["phylop-flank"]
data.head()

X = data.loc[data["annotation"] != "exonic", "constrained-score"].values.reshape(-1, 1)
print(X.shape)
print(X)
# fit 10 2 componenent gaussian mixture models
models = [None for i in range(10)]
for i in range(10):
    means_init = np.array([[0],[1]])
    weights_init = np.array([.9,.1])
    models[i] = GaussianMixture(2,
                                covariance_type="spherical",
                                random_state=i+1,
                                means_init=means_init, weights_init=weights_init).fit(X)
    
# compute the AIC and the BIC
AIC = [m.aic(X) for m in models]
BIC = [m.bic(X) for m in models]

fig, ax = plt.subplots(2,2, tight_layout=True, figsize=(12,8))
M_best = models[np.argmin(BIC)]
m1 = M_best.means_[0][0]
m2 = M_best.means_[1][0]
std1 = M_best.covariances_[0]
std2 = M_best.covariances_[1]
x = np.linspace(-10, 10, 1000)
logprob = M_best.score_samples(x.reshape(-1, 1))
responsibilities = M_best.predict_proba(x.reshape(-1, 1))
pdf = np.exp(logprob)
pdf_individual = responsibilities * pdf[:, np.newaxis]

ax[0,0].hist(X, 100, density=True, histtype='stepfilled', alpha=0.4)
ax[0,0].plot(x, pdf, '-k')
ax[0,0].plot(x, pdf_individual, '--k')
ax[0,0].set_xlim([-10,10])
ax[0,0].text(0.04, 0.96, "Best model",
        ha='left', va='top', transform=ax[0,0].transAxes)

ax[0,0].text(0.96, 0.96, "Parameters\nmu 1 = {0:.2f}\nmu 2 = {1:.2f}\nstd 1 = {2:.2f}\nstd 2 = {3:.2f}".format(m1,m2, std1, std2),
        ha='right', va='top', transform=ax[0,0].transAxes)

y_min, y_max = ax[0,0].get_ylim()
ax[0,0].vlines([m1,m2], y_min, y_max, linestyles="dashed", color="red")

# ax[0,0].set_xlabel('Constrained score')
ax[0,0].set_xlabel('Core PhyloP')
ax[0,0].set_ylabel('$p(x)$')
ax[0,0].set_title(tf)
x0,x1 = ax[0,0].get_xlim()
y0,y1 = ax[0,0].get_ylim()
ax[0,0].set_aspect(abs(x1-x0)/abs(y1-y0))

ax[0,1].hist(X, 100, density=True, histtype='stepfilled', alpha=0.4)
for color_i, i  in enumerate(np.argsort(BIC)[::-1]):
    M = models[i]
    x = np.linspace(-10, 10, 1000)
    logprob = M.score_samples(x.reshape(-1, 1))
    responsibilities = M.predict_proba(x.reshape(-1, 1))
    pdf = np.exp(logprob)
#     pdf_individual = responsibilities * pdf[:, np.newaxis]
    if color_i == 9:
        ax[0,1].plot(x, pdf, color="k", label = "rank: {}, BIC: {}".format(10-color_i, round(BIC[i],0)))
    else:
        ax[0,1].plot(x, pdf, color=cm.tab10(color_i/10), label = "rank: {}, BIC: {}".format(10-color_i, round(BIC[i],0)))
#     ax[0].plot(x, pdf_individual, '--k')

ax[0,1].set_xlim([-10,10])
ax[0,1].set_xlabel('Core PhyloP')
ax[0,1].set_ylabel('$p(x)$')
ax[0,1].set_title(tf)
ax[0,1].legend()
    
x = np.array([.1*x for x in range(-100,100)])
probs = M_best.predict_proba(x.reshape(-1, 1))
if m1 < m2:
    ax[1,0].plot(x, probs[:,0], "-k", label="m1")
    ax[1,0].plot(x, probs[:,1], "--k", label="m2")
else:
    ax[1,0].plot(x, probs[:,1], "-k", label="m1")
    ax[1,0].plot(x, probs[:,0], "--k", label="m2")
ax[1,0].set_xlabel("Core PhyloP")
ax[1,0].set_ylabel("Probability")
ax[1,0].legend()

Ps = [0.5,0.6,0.7,0.8,0.9,0.95]
n_constrained = []
thresholds = []
for index, p in enumerate(Ps):
    if m1 < m2:
        idx1 = np.argwhere(probs[:,1] > p).flatten().tolist()
    else:
        if std2 < std1:
            idx1 = np.argwhere(probs[:,0] > p).flatten().tolist()
        else:
            idx1 = np.argwhere(probs[:,1] > p).flatten().tolist()
    
    idx2 = np.argwhere(x > 0).flatten().tolist()
    for i in idx2:
        if i in idx1:
            break

    thresh = x[i]
    thresholds.append(thresh)
    n_constrained.append(len([x for x in X if x > thresh]))
    ax[1,1].bar(index+1, len([x for x in X if x > thresh]),edgecolor='black', color='None')
ax[1,1].set_xticks([i+1 for i in range(len(Ps))])
ax[1,1].set_xticklabels(Ps)
ax[1,1].text(0.96, 0.96, "N total = {}".format(X.shape[0]),ha='right', va='top', transform=ax[1,1].transAxes)
ax[1,1].set_xlabel("Probability threshold")
ax[1,1].set_ylabel("# TFBSs")
plt.savefig("./GMM-Core-Only/" + tf + ".GMM.pdf")

with open("./GMM-Core-Only/" + tf + ".txt", "w") as f:
    if m1 < m2:
        print(tf, m1, m2, std1, std2, BIC[np.argmin(BIC)], *n_constrained, n_p, n_d, file=f, sep="\t")
    else:
        print(tf, m2, m1, std2, std1, BIC[np.argmin(BIC)], *n_constrained, n_p, n_d, file=f, sep="\t")

print("Threshold = {}".format(thresholds[0]))

data["constrained"] = False
data.loc[(data["constrained-score"] > thresholds[0]), "constrained"] = True

probs = M_best.predict_proba(data["constrained-score"].values.reshape(-1, 1))
print(probs.shape)
print(probs[0,:])
data["p"] = probs[:,1]
data.loc[data["constrained-score"] < 0, "p"] = 0
# data.to_csv("./7-Constrained-TFBSs-Core-Only/" + tf + ".bed",
#             sep="\t",
#             header=False,
#             index=False,
#             float_format='%.15f')

data.to_csv("./tmp2/" + tf + ".bed",
            sep="\t",
            header=False,
            index=False,
            float_format='%.5f')

fig, ax = plt.subplots(1,2, figsize=(16,6), 
                       tight_layout=True, gridspec_kw={'width_ratios': [1, 2]})

M_best = models[np.argmin(BIC)]
m1 = M_best.means_[0][0]
m2 = M_best.means_[1][0]
std1 = M_best.covariances_[0]
std2 = M_best.covariances_[1]
x = np.linspace(-10, 10, 1000)
logprob = M_best.score_samples(x.reshape(-1, 1))
responsibilities = M_best.predict_proba(x.reshape(-1, 1))
pdf = np.exp(logprob)
pdf_individual = responsibilities * pdf[:, np.newaxis]

ax[0].hist(X, 100, density=True, histtype='stepfilled', alpha=0.4)
ax[0].plot(x, pdf, '-k')
ax[0].plot(x, pdf_individual, '--k')
ax[0].set_xlim([-10,10])
ax[0].text(0.04, 0.96, "Best model",
        ha='left', va='top', transform=ax[0].transAxes)

ax[0].text(0.96, 0.96, "Parameters\nmu 1 = {0:.2f}\nmu 2 = {1:.2f}\nstd 1 = {2:.2f}\nstd 2 = {3:.2f}".format(m1,m2, std1, std2),
        ha='right', va='top', transform=ax[0].transAxes)

y_min, y_max = ax[0].get_ylim()
ax[0].vlines([m1,m2], y_min, y_max, linestyles="dashed", color="red")

# ax[0,0].set_xlabel('Constrained score')
ax[0].set_xlabel('Core PhyloP')
ax[0].set_ylabel('$p(x)$')
ax[0].set_title(tf)

# x0,x1 = ax[0].get_xlim()
# y0,y1 = ax[0].get_ylim()
# ax[0].set_aspect(abs(x1-x0)/abs(y1-y0))

# Hide the right and top spines
ax[0].spines['right'].set_visible(False)
ax[0].spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
ax[0].yaxis.set_ticks_position('left')
ax[0].xaxis.set_ticks_position('bottom')


pfm = np.sum(data["seq"].map(encode_sequence))
ppm = pfm / np.sum(pfm, axis=1, keepdims=True)

indices = data.index[(data["annotation"] != "exonic") & (data["constrained"] == True) & (data["strand"] == "+")]
y1 = np.mean(phylop[indices], axis=0)
indices = data.index[(data["annotation"] != "exonic") & (data["constrained"] == True) & (data["strand"] == "-")]
y2 = np.mean(phylop[indices], axis=0)[::-1]
constrained = np.mean((y1,y2), axis=0)

indices = data.index[(data["annotation"] != "exonic") & (data["constrained"] == False) & (data["strand"] == "+")]
y1 = np.mean(phylop[indices], axis=0)
indices = data.index[(data["annotation"] != "exonic") & (data["constrained"] == False) & (data["strand"] == "-")]
y2 = np.mean(phylop[indices], axis=0)[::-1]
not_constrained = np.mean((y1,y2), axis=0)

x = np.arange(50+w)
ax[1].plot(x, constrained[75:125+w], color="red")
ax[1].plot(x, not_constrained[75:125+w], color="blue")
ax[1].set_xticks([0,25-.5,25+w-0.5,50+w])
ax[1].set_xticklabels(["-25bp","start","end","+25bp"])
bottom, top = ax[1].get_ylim() 
ax[1].set_ylim(0-.4*top, top)
ax[1].set_ylabel("241 way phylop")
inset_width = str(w / (50+w) * 100) + "%"
inset = inset_axes(ax[1],
                        width=inset_width, # width = 30% of parent_bbox
                        height=.5, # height : 1 inch
                        loc=8)

ppm = pd.DataFrame(ppm, columns=["A", "C", "G", "T"])
ppm = ppm.set_index(ppm.index + 1)

inset_logo = logomaker.Logo(ppm.applymap(get_information_content),
                           ax=inset)
inset.axis('off')

# Hide the right and top spines
ax[1].spines['right'].set_visible(False)
ax[1].spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
ax[1].yaxis.set_ticks_position('left')
ax[1].xaxis.set_ticks_position('bottom')

plt.savefig("./GMM-Core-Only/" + tf + ".phylop.pdf")