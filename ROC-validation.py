import numpy as np
from sklearn import metrics as skm
import argparse
import random
import math
import os
import csv
import matplotlib.pyplot as plt

# Setting a random seed at the start of a program can
# be a good habit to get into, because it removes one
# obstacle to reproducible data.
random.seed(10000)

def parse_arguments():
    # Here is where the input arguments are parsed.
    # There are two arguments: -e, which is the 
    # folder in which program is being run, and -f,
    # which corresponds to the --fdrlimit argument to
    # Wellington (and is set arbitrarily, often to -5,
    # for pipelines invoking HINT).
    
    parser = argparse.ArgumentParser(description='Produces a ROC plot from ChIP-seq validation data.')
    parser.add_argument('-e','--experiment',type=str,help='the output directory')
    parser.add_argument('-f','--fdrlimit',type=int,help='the log10 p-value cutoff')
    
    return parser.parse_args()


def read_tfpn_bedfile(file):
    
    # This is a function that reads a .bed file
    # of the PWM hits AKA events and relays
    # that to a list object.
    
    # Here we open the .bed file.
    tfpn_bedfile = open(file,'rU')
    
    # Here we initialize the returned list.
    tfpn_list = []
    
    # Here we iterate through it, pulling each tab-
    # delimited line as a list into our list.
    for line in csv.reader(tfpn_bedfile, delimiter='\n'):
        for coordinate in csv.reader(line, delimiter='\t'):
            tfpn_list.append(coordinate)

    # Closing files is hygienic.
    tfpn_bedfile.close()

    return tfpn_list


def find_file_length(file):
    
    # This is a function that reads a .bed file
    # of the open chromatin intervals and returns
    # the number of lines in a file.
    
    # Here we open the .bed file.
    file_to_check = open(file,'rU')
    
    # Here we read the line number into a length.
    length = len(file_to_check.readlines())

    # Closing files is hygienic.
    file_to_check.close()
    
    return length


def generate_lognormal(limit):

    # This is a very simple function that takes a
    # lower limit and generates a lognormal value
    # above that limit.

    # Here we start the process with a lognormal
    # value.
    lognorm_value = np.random.lognormal()
    
    # We repeat this process until we have a value
    # above our cutoff.
    while (math.pow(10.0,-lognorm_value) < math.pow(10.0,limit)):
        lognorm_value = np.random.lognormal()

    # This process functions like sampling from a
    # lognormal distribution that's been sliced at
    # a certain point.
    return math.pow(10.0,-lognorm_value)


def produce_pred_values(directory, limit):

    # You don't care about the events as much as
    # how many positives and negatives there are.
    tp_len = find_file_length(directory+'/true_positives.bed')
    fp_len = find_file_length(directory+'/false_positives.bed')
    fn_len = find_file_length(directory+'/false_negatives.bed')
    tn_len = find_file_length(directory+'/true_negatives.bed')

    # For the true calls, there's a p-value given
    # by the footprinting algorithm that it helps
    # to hold onto, so you read those in as lists.
    tp_list = read_tfpn_bedfile(directory+'/true_positives.bed')
    fp_list = read_tfpn_bedfile(directory+'/false_positives.bed')

    # Initializing numpy arrays for each event type.
    tp_values = np.zeros(tp_len)
    fp_values = np.zeros(fp_len)
    fn_values = np.zeros(fn_len)
    tn_values = np.zeros(tn_len)

    # High scores go to low p-values. The actual
    # number isn't important for a ROC plot, just
    # the order. The first guesses have to have the
    # highest values associated with them.
    for i in range(len(tp_values)):
        tp_values[i] = 1.0 - math.pow(10.0,-1*abs(float(tp_list[i][14])))
    for i in range(len(fp_values)):
        fp_values[i] = 1.0 - math.pow(10.0,-1*abs(float(fp_list[i][14])))

    # The falses, lacking a value, are given random
    # values.
    for i in range(len(fn_values)):
        fn_values[i] = generate_lognormal(limit)
    for i in range(len(tn_values)):
        tn_values[i] = generate_lognormal(limit)

    return tp_values,fp_values,fn_values,tn_values


def craft_predlabel_vectors(tpv,fpv,fnv,tnv):

    # Initializing numpy arrays for all events.
    predictions = np.zeros(len(tpv)+len(tnv)+len(fpv)+len(fnv))
    labels = np.zeros(len(tpv)+len(tnv)+len(fpv)+len(fnv))
    
    # The predictions and labels arrays are
    # edited to correspond with individual events
    # and their true/false positivity/negativity.
    for j in range(len(tpv)):
        predictions[j] = tpv[j]
        labels[j] = 1
    for j in range(len(fpv)):
        predictions[j+len(tpv)] = fpv[j]
        labels[j+len(tpv)] = 0
    for j in range(len(fnv)):
        predictions[j+len(tpv)+len(fpv)] = fnv[j]
        labels[j+len(tpv)+len(fpv)] = 1
    for j in range(len(tnv)):
        predictions[j+len(tpv)+len(fpv)+len(fnv)] = tnv[j]
        labels[j+len(tpv)+len(fpv)+len(fnv)] = 0

    return predictions, labels


def iterate_through_chip(experiment, limit):
    
    # Initializing fpr vector.
    mean_tpr = []
    mean_fpr = np.linspace(0,1,100)
    j = 0
    
    # Iterates through the folder that contains
    # the event bedfiles.
    for chip_target in os.listdir(experiment+'/ChIP-validation'):
        # Produces the lists/lengths for ROC analysis.
        tpv,fpv,fnv,tnv = produce_pred_values(experiment
                                              +'/ChIP-validation/'
                                              +chip_target, limit)
        # Checks to make sure the lists aren't all
        # true negatives.
        if len(tpv)+len(fpv)+len(fnv) > 0:
            # Produces the vectors from the event lists.
            preds, labels = craft_predlabel_vectors(tpv,fpv,fnv,tnv)
            # Performs the ROC analysis.
            fpr,tpr,thresholds = skm.roc_curve(labels,preds)
            new_mean = np.interp(mean_fpr,fpr,tpr).tolist()
            if j == 0: mean_tpr += new_mean
            else:
                for point in range(len(mean_tpr)):
                    mean_tpr[point] += new_mean[point]
            mean_tpr[0] = 0.0
            roc_auc = skm.auc(fpr, tpr)
            print chip_target,"AUC:",roc_auc
            j += 1

    for point in range(len(mean_tpr)):
        mean_tpr[point] /= j
    mean_tpr[-1] = 1.0
    mean_auc = skm.auc(mean_fpr,mean_tpr)
    print "Mean AUC:",mean_auc


def __main__():

    # This function governs all other functions and runs
    # only when the program is run directly.

    # Here the input arguments are accepted. These arguments
    # are passed to variables so that the entire args object
    # doesn't have to be passed to each function.
    args = parse_arguments()
    experiment = args.experiment
    limit = args.fdrlimit

    # Here we read the input and produce figures.
    iterate_through_chip(experiment,limit)


# The condition that the program is run directly.
# This makes unit testing easy.
if __name__ == '__main__':
    __main__()