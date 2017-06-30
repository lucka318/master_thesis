from Bio import SeqIO
import readerrors as re

fasta_sequences = SeqIO.parse(open("../data/start_data/scerevisiae_reference.fasta"),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.description, str(fasta.seq)
	print(name)
	#print(fasta)
#re.readerrors("data/start_data/scerevisiae_reference.fasta", "data/sam_files/ecoli_pb_1.sam")