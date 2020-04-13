import cPickle as pic

seq_unique = pic.load(open('seq_unique.pkl','rb'))
index_unique = pic.load(open('index_unique.pkl','rb'))
seq_record_dict = pic.load(open('seq_record_dict.pkl','rb'))
index_record_dict = pic.load(open('index_record_dict.pkl','rb'))
target_seq = 'CGAGTAAT'

def hamming(seq1,seq2):
	return sum(base1 != base2 for base1,base2 in zip(seq1,seq2))
index_hamming=[]
for ID in seq_unique:
	index = seq_record_dict[ID].description.split(':')[-1]
	sequence =  seq_record_dict[ID].seq
	index_hamming.append(hamming(index,target_seq))
	print index, target_seq, hamming(index,target_seq)
seq_hamming=[]
for ID in index_unique:
	index = seq_record_dict[ID].description.split(':')[-1]
	sequence =  seq_record_dict[ID].seq
	seq = sequence[-9:-1]
	seq_hamming.append(hamming(seq,target_seq))
	print index, target_seq, hamming(seq,target_seq)


plt.hist(index_hamming)
plt.show()
plt.hist(seq_hamming)
plt.show()
