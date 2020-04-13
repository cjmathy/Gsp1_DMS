import sys
import cPickle as pic 
import numpy as np 
import matplotlib.pyplot as plt 

matrix = pic.load(open(sys.argv[1],'rb'))
matrix = np.ma.masked_invalid(matrix)
transform = {'A': 9, 'C': 8, 'E': 20, 'D': 19, 'G': 10, 'F': 2, 'I': 5, 'H': 16, 'K': 18, '*': 0, 'M': 6, 'L': 4, 'N': 14, 'Q': 15, 'P': 11, 'S': 12, 'R': 17, 'T': 13, 'W': 1, 'V': 7, 'Y': 3}
WT_seq = 'MSAPAANGEVPTFKLVLVGDGGTGKTTFVKRHLTGEFEKKYIATIGVEVHPLSFYTNFGEIKFDVWDTAGQEKFGGLRDGYYINAQCAIIMFDVTSRITYKNVPNWHRDLVRVCENIPIVLCGNKVDVKERKVKAKTITFHRKKNLQYYDISAKSNYNFEKPFLWLARKLAGNPQLEFVASPALAPPEVQVDEQLMQQYQQEMEQATALPLPDEDDADL*'



def stops(matrix):
	return matrix[0,]

def syn(matrix):
	WT_syn_scores = []
	for i, letter in enumerate(WT_seq):
		position = i
		mutant_number= int(transform[letter])
		score = matrix[mutant_number][position]
		WT_syn_scores.append(score)
	return WT_syn_scores

WT_syn = np.asarray(syn(matrix))
stop_scores = stops(matrix)
print WT_syn
WT_syn_mean =  0

plt.figure()
plt.suptitle("Fitness Score Distribution")
bins =  np.arange(np.nanmin(matrix), np.nanmax(matrix), 0.1)
plt.hist((matrix[~np.isnan(matrix)]-WT_syn_mean), bins = bins, alpha = 0.5, color = "b", label = "All")
plt.hist(WT_syn[~np.isnan(WT_syn)]-WT_syn_mean, bins = bins, alpha = 0.5, color = "g", label = "Synonymous")
plt.hist(stop_scores[~np.isnan(stop_scores)]-WT_syn_mean, bins = bins, alpha = 0.75, color = "r", label = "Nonsense")

plt.xlabel("Fitness Score")
plt.ylabel("Count")
plt.legend()

### Plotting score vs position for nonsense mut
# plt.figure()
# plt.suptitle("Fitness Score by Position")
# plt.plot(syn_masked, "go", ms = 5, label = "Synonymous")
# plt.plot(stop,"ro", ms = 5, label = "Nonsense")

# plt.xlabel("Position")
# plt.ylabel("Fitness Score")
# plt.legend()

### Show plot
plt.savefig('%s_hist.png'%(sys.argv[1][:-4]), dpi =300)
#plt.show()