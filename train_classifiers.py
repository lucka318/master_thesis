import argparse
import time
import sys
import numpy as np
import readerrors as re
import ml


def main(arguments):

	parser = argparse.ArgumentParser(description=__doc__)

	parser.add_argument("-t", "--train",
		required="True",
		dest="train",
		help="CSV file.")
	parser.add_argument("-n", "--number",
		required="True",
		dest="num",
		help="Number of training examples.")

	args = parser.parse_args(arguments)

	train_file = open(args.train)
  	train_data = np.genfromtxt(train_file, dtype=int, delimiter=',')
	ml.trainSVM(train_data,int(args.num))
	ml.trainRFC(train_data,int(args.num))
	ml.trainGBC(train_data,int(args.num))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))