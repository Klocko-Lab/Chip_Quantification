# dataconvert.py: Convert from  homer file format to N by N matrix
# OluwadareLab, 2021
# Usage information: python dataconvert.py -h

import numpy as np
import argparse

from array import array

def main(data_loc, seq_file):

	with open(data_loc, 'r') as fl:
		col_data = fl.readline().split("\t")[:-1]
		n_cols = len(col_data)

	data = np.genfromtxt(data_loc, delimiter="\t", dtype = None, skip_header=1, usecols=np.arange(2, n_cols))

	file = open(seq_file, 'w')

	for seq in data:
		for value in seq:
			file.write(str(value))
			file.write(" ")
		file.write("\n")
	file.close()


if __name__ == '__main__':
	parser=argparse.ArgumentParser()
	parser.add_argument("hicExplorer_out",help="homer file generated from hicExplorer")
	parser.add_argument("contact_matrix",help="File name of a output white space delimited square contact matrix")

	args=parser.parse_args()
	if args.hicExplorer_out:
		data_loc = args.hicExplorer_out # this is the  v (h5 file is converted to homer
	if args.contact_matrix:
		seq_file = args.contact_matrix  # this is the n x n matrix output in a txt file
		
	main(data_loc, seq_file)