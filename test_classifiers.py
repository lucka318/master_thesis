import argparse
import time
import sys
import numpy as np
import readerrors as re
import ml


def main(arguments):

	parser = argparse.ArgumentParser(description=__doc__)

	parser.add_argument("-t", "--test",
		required="True",
		dest="test",
		help="CSV file.")
	parser.add_argument("-n", "--number",
		required="True",
		dest="num",
		help="Number of test examples.")

	args = parser.parse_args(arguments)

	test_file = open(args.test)
  	test_data = np.genfromtxt(test_file, dtype=int, delimiter=',')
	ml.testSVM(test_data,int(args.num))
	ml.testRFC(test_data,int(args.num))
	ml.testGBC(test_data,int(args.num))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))