import csv
import sys
csv.field_size_limit(sys.maxsize)

import glob
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-f','--folder',type=str,help='the folder')
args = parser.parse_args()

folder = args.folder

file = open(folder+'/peaks.bedgraph')
bedgraph = []
for line in csv.reader(file,delimiter='\n'):
	for tab in csv.reader(line,delimiter='\t'):
		bedgraph.append(tab)

peak_heights = []
for line in bedgraph:
	height = float(line[3])
	peak_heights.append(height)

slurms = glob.glob('slurm-*')

mapped_reads = {}
for slurm in slurms:
	openfile = open(slurm).read()
	readsend = openfile.find(' reads; of these:')
	folderend = openfile.find('/picard.bam as a bam file')
	if readsend != -1:
		key = openfile[:folderend].split('/')[-1]
		mapped_reads[key] = int(openfile[:readsend].split('\n')[-1])
		print key,int(openfile[:readsend].split('\n')[-1])

return_bedgraph = []
for line in bedgraph:
	line = line[:3] + [1000000*float(line[3])/mapped_reads[folder]]
	return_bedgraph.append(line)

file = open('return/'+folder+'.bedgraph','w')
for line in return_bedgraph:
	file.write(line[0]+'\t'+line[1]+'\t'+line[2]+'\t'+str(line[3])+'\n')