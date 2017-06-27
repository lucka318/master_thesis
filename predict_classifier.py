import argparse
import time
import sys
import numpy as np
import readerrors as re
import ml


def main(arguments):

	parser = argparse.ArgumentParser(description=__doc__)

	parser.add_argument("-p", "--predict",
		required="True",
		dest="predict",
		help="CSV file.")

	args = parser.parse_args(arguments)
	predict_file = open(args.predict)
	predict_data = np.genfromtxt(predict_file, dtype=int, delimiter=',')
	predict_data = ml.predictSVM(predict_data)
	with open("predict_SVM.csv",'w') as f:
		np.savetxt(f, predict_data.astype(int), fmt='%i', delimiter=',')

	predict_data = ml.predictRFC(predict_data)
	with open("predict_RFC.csv",'w') as f:
		np.savetxt(f, predict_data.astype(int), fmt='%i', delimiter=',')


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))