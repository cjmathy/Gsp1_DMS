import sys
import numpy as np 
import cPickle as pic 
import matplotlib.pyplot as plt
from scipy.stats import norm
import cPickle as pic

matrix = pic.load(open(sys.argv[1]))
outstring = sys.argv[1][:-4]
bins = np.linspace(np.nanmin(matrix),np.nanmax(matrix), 101)
transform = {'A': 9, 'C': 8, 'E': 20, 'D': 19, 'G': 10, 'F': 2, 'I': 5, 'H': 16, 'K': 18, '*': 0, 'M': 6, 'L': 4, 'N': 14, 'Q': 15, 'P': 11, 'S': 12, 'R': 17, 'T': 13, 'W': 1, 'V': 7, 'Y': 3}
transform_inv = {0: 'STOP', 1: 'W', 2: 'F', 3: 'Y', 4: 'L', 5: 'I', 6: 'M', 7: 'V', 8: 'C', 9: 'A', 10: 'G', 11: 'P', 12: 'S', 13: 'T', 14: 'N', 15: 'Q', 16: 'H', 17: 'R', 18: 'K', 19: 'D', 20: 'E'}
p = 0.01
Z = norm.ppf(1-(p)/2)
WT_seq = 'MSAPAANGEVPTFKLVLVGDGGTGKTTFVKRHLTGEFEKKYIATIGVEVHPLSFYTNFGEIKFDVWDTAGQEKFGGLRDGYYINAQCAIIMFDVTSRITYKNVPNWHRDLVRVCENIPIVLCGNKVDVKERKVKAKTITFHRKKNLQYYDISAKSNYNFEKPFLWLARKLAGNPQLEFVASPALAPPEVQVDEQLMQQYQQEMEQATALPLPDEDDADL*'
stops = matrix[0, :][:-20]
stops = stops[np.isfinite(stops)]
stops_mu, stops_sigma = norm.fit(stops)
WT_syn = []
for i, aa in enumerate(WT_seq):
    pos = i
    aa_number = transform[aa]
    fitness = matrix[aa_number , pos]  
    WT_syn.append(fitness)
WT_syn = np.array(WT_syn)
WT_syn = WT_syn[np.isfinite(WT_syn)]
syn_mu, syn_sigma = norm.fit(WT_syn)
Zsyn = (matrix - syn_mu) / syn_sigma
Zstops = (matrix - stops_mu) / stops_sigma
WT_like = matrix[np.where((Zsyn > -Z) & (Zsyn < Z))]
matrix[np.where((Zsyn > -Z) & (Zsyn < Z))] = 0.0
better = matrix[np.where((Zsyn >= Z))]
matrix[np.where((Zsyn >= Z))] = 0.5
stop_like = matrix[np.where((Zstops < Z))]
matrix[np.where((Zstops < Z))] = -1.0
inter = matrix[np.where((Zstops >= Z) & (Zsyn <= - Z))]
matrix[np.where((Zstops >= Z) & (Zsyn <= - Z))] = -0.5


pic.dump(matrix, open('%s_bins.pkl'%(outstring),'wb'))

plot_s  = stop_like[np.isfinite(stop_like)].flatten()
plot_wt = WT_like[np.isfinite(WT_like)].flatten()
plot_a = better[np.isfinite(better)].flatten()
plot_i = inter[np.isfinite(inter)].flatten()

plt.figure()
# S = plt.hist(stop_like[np.isfinite(stop_like)].flatten(), bins = bins,color = 'r', alpha = 0.75, label = 'Stop Like')
# WT = plt.hist(WT_like[np.isfinite(WT_like)].flatten(), bins = bins,color = 'g', alpha = 0.75, label = 'WT Like')
# I = plt.hist(inter[np.isfinite(inter)].flatten(), bins = bins,color = 'orange', alpha = 0.75, label = 'Intermediate')
# A = plt.hist(better[np.isfinite(better)].flatten(), bins = bins,color = 'blue', alpha = 0.75, label = 'Advantageous')
stacked = plt.hist([plot_s,plot_wt, plot_i, plot_a], normed = 0,bins = bins, color = ['r','g','orange','b'], label = ['Stop Like', 'WT Like', 'Intermediate', 'Advantageous'], stacked = 1)
plt.legend(loc = 0)
plt.xlabel('Fitness Score')
plt.ylabel('Counts')
#plt.savefig('%s_bins.png'%(outstring), dpi = 300)
plt.show()
