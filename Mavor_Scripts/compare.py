import sys
import cPickle as pic 
import numpy as np 
import matplotlib.pyplot as plt 

matrix1 = pic.load(open(sys.argv[1],'rb'))
matrix2 = pic.load(open(sys.argv[2],'rb'))
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

plt.figure()
plt.plot(matrix1.flatten(), matrix2.flatten(), 'k.', alpha = 0.75, mew = 0, ms = 10)
plt.plot(stops(matrix1), stops(matrix2), 'r.', alpha = 0.75, mew = 0, ms = 10)
plt.plot(syn(matrix1), syn(matrix2), 'g.', alpha = 0.75, mew = 0, ms = 10)
plt.xlabel(sys.argv[1][:-4])
plt.ylabel(sys.argv[2][:-4])
plt.savefig('%s_%s_scatter.png'%(sys.argv[1][:-4],sys.argv[2][:-4]), dpi = 300)
#plt.show()