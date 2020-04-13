import sys
import numpy as np 
import cPickle as pic

counts = open(sys.argv[1],'r')
counts.readline()
file_string = sys.argv[2]
input_counts_mat = np.empty([21,220])
sorted_counts_mat = np.empty([21,220])
input_counts_mat[:] = np.nan
sorted_counts_mat[:] = np.nan
transform = {'A': 9, 'C': 8, 'E': 20, 'D': 19, 'G': 10, 'F': 2, 'I': 5, 'H': 16, 'K': 18, '*': 0, 'M': 6, 'L': 4, 'N': 14, 'Q': 15, 'P': 11, 'S': 12, 'R': 17, 'T': 13, 'W': 1, 'V': 7, 'Y': 3}
WT_seq = 'MSAPAANGEVPTFKLVLVGDGGTGKTTFVKRHLTGEFEKKYIATIGVEVHPLSFYTNFGEIKFDVWDTAGQEKFGGLRDGYYINAQCAIIMFDVTSRITYKNVPNWHRDLVRVCENIPIVLCGNKVDVKERKVKAKTITFHRKKNLQYYDISAKSNYNFEKPFLWLARKLAGNPQLEFVASPALAPPEVQVDEQLMQQYQQEMEQATALPLPDEDDADL*'


def syn(matrix):
	WT_syn_scores = []
	for i, letter in enumerate(WT_seq):
		position = i
		mutant_number= int(transform[letter])
		score = matrix[mutant_number][position]
		WT_syn_scores.append(score)
	return WT_syn_scores

count_dict ={}
input_sum = 0
sorted_sum = 0

for line in counts:
	l_list = line.strip().split()
	allele = l_list[0]

 	counts = l_list[1:]
	try:
 		counts = [float(x) for x in counts]
 	except ValueError:
 		try:
 			counts = [float(counts[0]), 0.0]
 		except ValueError:
 			counts = [0.0,0.0]
 			allele = None
	input_sum+=counts[0]
	sorted_sum+=counts[1]
	try:
		pos = int(allele.split('_')[0])
		mut = allele.split('_')[1]
		input_counts_mat[transform[mut]][pos-1] = counts[0]
		sorted_counts_mat[transform[mut]][pos-1] = counts[1]
	except AttributeError:
		pass

input_sum = input_sum+0.5
sorted_sum = sorted_sum+0.5

input_counts_mat = input_counts_mat + 0.5
sorted_counts_mat = sorted_counts_mat + 0.5

freq_input_mat = input_counts_mat/input_sum
freq_sorted_mat = sorted_counts_mat/sorted_sum

score = np.log2(freq_sorted_mat) - np.log2(freq_input_mat)

WT_syn = syn(score)
WT_mean =  np.nanmean(WT_syn)

pic.dump(score - WT_mean , open('%s'%(file_string),'wb'))

