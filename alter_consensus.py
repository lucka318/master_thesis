import sys
import numpy as np
import pprint
import argparse
import csv
import numpy as np
from collections import defaultdict
from Bio import SeqIO

def genome_preprocessing(reference_file):

	with open(reference_file) as f:
		content = f.readlines()
	content = [x.strip() for x in content]
	genome=''.join(content[1:])
	return content[0], genome

def main(arguments):

	parser = argparse.ArgumentParser(description=__doc__)

	parser.add_argument("-c", "--consensus",
		required="True",
		dest="con",
		help="File containign consensus.")

	parser.add_argument("-s", "--csv",
    	required="True",
    	dest="csv",
    	help="CSV file with predicted values")

	parser.add_argument("-o", "--output",
    	required="True",
    	dest="out",
    	help="Fasta file for altered consensus.")

	args = parser.parse_args(arguments)

	f = open(args.csv)

	data = np.genfromtxt(f, dtype=int, delimiter=',')

	x_coord = data.shape[0]
	y_coord = data.shape[1]
	Y_ = data[:,y_coord-3:y_coord-2]
	homopolymer_coord = data[:,y_coord-2:y_coord]
	homopolymer_coord_dict = {}
	for i in range(0,x_coord):
		homopolymer_coord_dict[homopolymer_coord[i][0]] = Y_[i][0]

	fasta_sequences = SeqIO.parse(open(args.con),'fasta')

	l = open(args.out, 'w')
	indent = 0
	for seq in fasta_sequences:
		sequence = str(seq.seq)
		new_gen = ''
		letter = ''
		i = 0
		while(i < len(sequence)):
			letter = sequence[i]
			if(homopolymer_coord_dict.has_key(i + indent)):
				val = homopolymer_coord_dict[i + indent]
				for j in range(0,val):
					new_gen+=letter
				while(i < len(sequence) && sequence[i] == letter):
					i+=1
			else:
				new_gen+=sequence[i]
				i+=1
		indent += len(sequence)
		l.write('>' + seq.description + '\n')
		l.write(new_gen + '\n')

	l.close()
	print(len(new_gen))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
