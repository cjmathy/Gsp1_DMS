import sys
import cPickle as pic
import numpy as np

transform = {'A': 9, 'C': 8, 'E': 20, 'D': 19, 'G': 10, 'F': 2, 'I': 5, 'H': 16, 'K': 18, '*': 0, 'M': 6, 'L': 4, 'N': 14, 'Q': 15, 'P': 11, 'S': 12, 'R': 17, 'T': 13, 'W': 1, 'V': 7, 'Y': 3}

tsv =  open(sys.argv[1], 'r')
matrix_score = np.empty([21,220])
matrix_score[:] = np.nan
matrix_SE = np.empty([21,220])
matrix_SE[:] = np.nan
headers = tsv.readline().strip().split()
print headers
for line in tsv:
	l_list = line.strip().split()
	ID = l_list[0]
	score = float(l_list[1])
	SE = float(l_list[2])
	try:
		pos = int(ID.split('_')[0])
		mut = ID.split('_')[1]
		matrix_score[transform[mut]][pos-1] = score
		matrix_SE[transform[mut]][pos-1] = SE
	except ValueError:
		pass
pic.dump(matrix_score,open('%s_scores.pkl'%(sys.argv[2]),'w'))
pic.dump(matrix_SE,open('%s_SE.pkl'%(sys.argv[2]),'w'))