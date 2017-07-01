import argparse
import time
import sys
import numpy as np
import readerrors as re


out_cons_set = "consensus_predict.csv"
out_freq_set = "consensus_error.csv"

def main(arguments):

	parser = argparse.ArgumentParser(description=__doc__)

	parser.add_argument("-r1", "--reference1",
		required="True",
		dest="ref1",
		help="File containign reference, fasta format.")

	# parser.add_argument("-r2", "--reference2",
	# 	required="True",
	# 	dest="ref2",
	# 	help="File containign reference, fasta format.")

	parser.add_argument("-s1", "--reads1",
    	required="True",
    	dest="reads1",
    	help="Reads aligned to the reference. SAM format.")

	# parser.add_argument("-s2", "--reads2",
 #    	required="True",
 #    	dest="reads2",
 #    	help="Reads aligned to the reference. SAM format.")

	parser.add_argument("-c", "--consensus",
		required="True",
		dest="cons",
		help="File containing consensus, fasta format.")

	parser.add_argument("-d", "--con_reads",
    	required="True",
    	dest="cons_reads",
    	help="Reads aligned to the consensus. SAM format.")

	args = parser.parse_args(arguments)

	read_max_1, ref_max_1, regions_dict_1, freq_dict_1 = re.readerrors(args.ref1, args.reads1)

	#read_max_2, regions_dict_2 = re.readerrors(args.ref2, args.reads2)

	read_max_c, ref_max_c, regions_dict_c, freq_dict_c = re.readerrors(args.cons, args.cons_reads)

	read_max = max(read_max_1, read_max_c, ref_max_1, ref_max_c)
	re.make_test_csv(regions_dict_1, read_max, "train.csv")
	#re.make_test_csv(regions_dict_2, read_max, "test.csv")
	re.make_test_csv(regions_dict_c, read_max, "predict.csv")

	re.make_freq_csv(freq_dict_1, read_max_1, ref_max_1, "freq_ref.csv")
	re.make_freq_csv(freq_dict_c, read_max_c, ref_max_c, "freq_cons.csv")


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
