import cPickle as pic
import sys

matrix = pic.load(open(sys.argv[1],'rb'))
outstring = sys.argv[1][:-4]
transform = {'A': 9, 'C': 8, 'E': 20, 'D': 19, 'G': 10, 'F': 2, 'I': 5, 'H': 16, 'K': 18, '*': 0, 'M': 6, 'L': 4, 'N': 14, 'Q': 15, 'P': 11, 'S': 12, 'R': 17, 'T': 13, 'W': 1, 'V': 7, 'Y': 3}
transform_inv = {0: '*', 1: 'W', 2: 'F', 3: 'Y', 4: 'L', 5: 'I', 6: 'M', 7: 'V', 8: 'C', 9: 'A', 10: 'G', 11: 'P', 12: 'S', 13: 'T', 14: 'N', 15: 'Q', 16: 'H', 17: 'R', 18: 'K', 19: 'D', 20: 'E'}


f = open('%s.csv'%(outstring),'w')
for i, resi in enumerate(matrix):
	mut = transform_inv[i]
	for pos, score in enumerate(resi):
		key = mut+str((pos+1))
		f.write('%s, %s\n' %(, score))
f.close
