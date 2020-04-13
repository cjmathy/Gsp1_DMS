import sys
from Bio import SeqIO as seq
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import scipy.spatial.distance as dist
import cPickle as pic

seq_reads = open(sys.argv[1],'r')
index_reads = open(sys.argv[2],'r')

target_seq = 'CGAGTAAT'

seq_record_dict = {}
index_record_dict ={}
# print record_dict
seq_miss = []
index_miss = []
seq_reads_id = []
index_reads_id = []
for record in seq.parse(seq_reads, 'fastq'):
	index = record.description.split(':')[-1]
	seq_record_dict[record.id] = record
	seq_reads_id.append(record.id)
# # 	if index == target_seq:
# 		pass
# 	else:
# 		if record.seq[-9:] == 'CGAGTAATA':
# 			seq_miss.append(record)

# print len(seq_miss)
for record in seq.parse(index_reads, 'fastq'):
	index = record.description.split(':')[-1]
	index_reads_id.append(record.id)
	index_record_dict[record.id] = record
	# if record.seq[-11:] == 'TACGAGTAATA':
	# 	pass
	# else:
	# 	if index == target_seq:
	# 		index_miss.append(record)
	# 	else:
# 	# 		print 'BAD'
seq_reads_id = set(seq_reads_id)
index_reads_id = set(index_reads_id)

seq_unique = seq_reads_id - index_reads_id
index_unique = index_reads_id - seq_reads_id 


index_hamming = []
seq_hamming = []
def hamming(seq1,seq2):
	return sum(base1 != base2 for base1,base2 in zip(seq1,seq2))

# for ID in seq_unique:
# 	index = seq_record_dict[ID].description.split(':')[-1]
# 	sequence =  seq_record_dict[ID].seq
# 	index_hamming.append(hamming(index,target_seq))
# 	print index, target_seq, hamming(index,target_seq)

check_me =[]

for ID in index_unique:
	#index = seq_record_dict[ID].description.split(':')[-1]
	sequence =  index_record_dict[ID].seq
	seq_string = sequence[-9:-1]
	seq_hamming.append(hamming(seq_string,target_seq))
	#print seq_string, target_seq, hamming(seq_string,target_seq)
	if hamming(seq_string,target_seq) > 2.0:
		check_me.append(ID)


# pic.dump(seq_unique, open('seq_unique.pkl','wb'))
# print 'one'
# pic.dump(index_unique, open('index_unique.pkl','wb'))
# print 'two'
# pic.dump(seq_record_dict, open('seq_record_dict.pkl','wb'))
# print 'three'
# pic.dump(index_record_dict, open('index_record_dict.pkl','wb'))

# print 'DONE'
# plt.figure()
# bins = [0,1,2,3,4,5,6,7,8,9]
# plt.suptitle('Index')
# plt.hist(index_hamming, bins = bins)
# plt.savefig('index_haming.png')
# plt.show()

# plt.suptitle('Sequence')
# plt.hist(seq_hamming, bins = bins)
# plt.savefig('seq_haming.png')
# plt.show()

# venn2([seq_reads_id, index_reads_id],['Grep on Read', 'Grep on Index'])
# plt.savefig('venn.png')

# plt.show()

indices = ['CGAGTAAT', 'TCTCCGGA','AATGAGCG', 'GGAATCTC',  'TTCTGAAT', 'ACGAATTC']

for ID in check_me:
	sequence =  index_record_dict[ID].seq
	seq_string = sequence[-9:-1]
	for index in indices:
		if hamming(seq_string,index) < 3.0:
			print seq_string, index, hamming(seq_string,index), ID

