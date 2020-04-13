import sys
import numpy as np 
import matplotlib.pyplot as plt 
import cPickle as pic 
from numpy import ma

parul =  pic.load(open(sys.argv[1],'rb'))
mavor = pic.load(open(sys.argv[2],'rb'))

#low_scores = np.where(parul <-7.5)
# print len(np.where(np.isnan(mavor[low_scores]))[0])

nan_scores = np.where(np.isnan(mavor))

# matrix = np.ma.masked_invalid(mavor[low_scores])

WT_seq = 'MSAPAANGEVPTFKLVLVGDGGTGKTTFVKRHLTGEFEKKYIATIGVEVHPLSFYTNFGEIKFDVWDTAGQEKFGGLRDGYYINAQCAIIMFDVTSRITYKNVPNWHRDLVRVCENIPIVLCGNKVDVKERKVKAKTITFHRKKNLQYYDISAKSNYNFEKPFLWLARKLAGNPQLEFVASPALAPPEVQVDEQLMQQYQQEMEQATALPLPDEDDADL*'
transform = {'A': 9, 'C': 8, 'E': 20, 'D': 19, 'G': 10, 'F': 2, 'I': 5, 'H': 16, 'K': 18, '*': 0, 'M': 6, 'L': 4, 'N': 14, 'Q': 15, 'P': 11, 'S': 12, 'R': 17, 'T': 13, 'W': 1, 'V': 7, 'Y': 3}
transform_inv = {0: 'STOP', 1: 'W', 2: 'F', 3: 'Y', 4: 'L', 5: 'I', 6: 'M', 7: 'V', 8: 'C', 9: 'A', 10: 'G', 11: 'P', 12: 'S', 13: 'T', 14: 'N', 15: 'Q', 16: 'H', 17: 'R', 18: 'K', 19: 'D', 20: 'E'}

matrix = parul.T

mut = np.where(np.logical_and(matrix <-7.5, np.isnan(mavor.T)))[1]
pos = np.where(np.logical_and(matrix <-7.5, np.isnan(mavor.T)))[0]
l = []
for i, foo in enumerate(pos):
	if mut[i] != 0:
		#print WT_seq[foo]+str(foo+1)+transform_inv[mut[i]]
		l.append(WT_seq[foo]+str(foo+1)+transform_inv[mut[i]])

print l



# plt.figure()
# plt.hist(matrix, bins = 75)
# plt.xlabel('Fitness Sscore')
# plt.ylabel('Counts')
# plt.savefig('Low_scores_Mavor_data.png', dpi = 300)
# plt.show()

