import csv
import sys
csv.field_size_limit(sys.maxsize)
import argparse
import pickle
import numpy as np


def parse_args():

	parser = argparse.ArgumentParser(description='make .npy files for deep learning')
	parser.add_argument('-f','--fasta',type=str,help='the fasta file, with headers')
	parser.add_argument('-o','--output',type=str,help='the measured output of the system, as a pickled dict')

	return parser.parse_args()


def read_in_fasta(inpath):

	openfile = open(inpath)
	returndict = {}
	fastaname = ''; fastaseq = ''
	for line in csv.reader(openfile,delimiter='\n'):
		if line[0][0] == '>':
			returndict[fastaname] = fastaseq
			fastaname = line[0][1:]
			fastaseq = ''
		else:
			fastaseq += line[0]

	return returndict


def read_in_pickle(inpath):

	with open(inpath,'rb') as openfile:
		returndict = pickle.load(openfile)
	return returndict


def convert_table(table):

	returndict = {}
	for line in table[1:]:
		ENST = line[0].split('.')[0]
		FPKM = float(line[6])
		returndict[ENST] = FPKM

	return returndict


def pickle_(object,outpath):

	with open(outpath,'wb') as openfile:
		pickle.dump(object,openfile)


def DNA_1hE(sequence_DNA):

	string = sequence_DNA.upper()
	_1hElol = []
	for letter in string:
		if letter == 'A':
			_1hElol.append([1,0,0,0])
		elif letter == 'C':
			_1hElol.append([0,1,0,0])
		elif letter == 'G':
			_1hElol.append([0,0,1,0])
		elif letter == 'T':
			_1hElol.append([0,0,0,1])
		else:
			_1hElol.append([0,0,0,0])

	return _1hElol


def link_together(fasta,ydict):

	x_list = []; y_list = []
	for key in ydict:
		if key in fasta:
			sequence = fasta[key]
			one_hot = DNA_1hE(sequence)
			output = ydict[key]
			x_list.append(output)
			y_list.append(one_hot)
	x_arr = np.array(x_list)
	y_arr = np.array(y_list)
	return x_arr,y_arr


def print_nparrays(x_array,y_array):

	np.save('x.npy',x_array)
	np.save('y.npy',y_array)


def __main__():
	
	args = parse_args()
	fastapath = args.fasta; outputpath = args.output
	fasta = read_in_fasta(fastapath)
	ydict = read_in_pickle(outputpath)
	x_arr,y_arr = link_together(fasta,ydict)
	print_nparrays(x_arr,y_arr)


if __name__ == '__main__':
	__main__()