import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import csv
import sys
csv.field_size_limit(sys.maxsize)

import numpy as np
import pandas as pd
import sklearn.naive_bayes as sknb
import sklearn.ensemble as ske
from math import pow


def read_ylist(path):
    
    # This is a function that reads a .txt file
    # of the chrY gene list and relays that
    # to a list object.

    ylist = open(path,'rU')
    
    return_list = []
    
    for line in csv.reader(ylist, delimiter='\n'):
        return_list += line

    # Closing files is hygienic.
    ylist.close()
    
    return return_list


def read_converter(path):
    
    # The concept here is to produce
    # a dictionary to convert values
    # from one regime to another.

    converter = open(path,'rU')
    
    return_dict = {}
    
    # Only read in lines that match a
    # certain format.
    for line in csv.reader(converter, delimiter='\n'):
        for tab in csv.reader(line,delimiter='\t'):
            if len(tab) > 13:
                return_dict[tab[13].lower()] = tab[12]
    
    # Closing files is hygienic.
    converter.close()
    
    return return_dict


def read_datafile(path):

    # Reads the data into a dictionary format
    # where each ENSG is the key to its
    # expression.
    
    data = open(path,'rU')
    
    return_dict = {}
    
    for line in csv.reader(data, delimiter='\n'):
        for tab in csv.reader(line,delimiter='\t'):
            return_dict[tab[0]] = [float(value) for value in tab[1:]]
    
    # Closing files is hygienic.
    data.close()
    
    return return_dict


def print_file(path,ylist,converter,data,ratio=False):

    # These are the lazily-coded values for each of 
    # these standards.

    printable_list = [['NAME','ID_REF','SAMPLE 5','SAMPLE 6','SAMPLE 7','SAMPLE 8','SAMPLE 1','SAMPLE 2','SAMPLE 3','SAMPLE 4','SAMPLE 9','SAMPLE 10','SAMPLE 11','SAMPLE 12','PSG101']]
    
    # GAPDH expression in each sample.
    gapdh = [pow(value,2) for value in [-0.15707731,-0.42139864,0.04047537,-0.34588385,0.24991941,-0.04047537,0.044377804,-0.12586355,0.34507227,0.2800212,0.20464182,-0.0790391,-1.2639866]]
    
    # ACTB expression in each sample.
    actb = [pow(value,2) for value in [-0.1605568,-0.07634735,1.016263,0.07634735,0.8745966,-0.24687767,-0.29980946,-0.29821014,-0.27721214,0.13205433,0.9464884,0.5554552,-0.12448788]]
    
    # Iterate through list of chrY genes...
    for ygene in ylist:
        # Check if the gene is in the dictionary...
        if ygene.lower() in converter:
            # Initialize a list and add both names
            # for the gene...
            list_to_append = [ygene,converter[ygene.lower()]]
            # Convert back from log2 space...
            values = [pow(value,2) for value in data[converter[ygene.lower()]]]
            # (Here's an optional step to rescale by
            # ACTB expression.)
            if ratio:
                values = [values[i]/actb[i] for i in range(len(actb))]
            list_to_append += values
            # And append the list!
            printable_list.append(list_to_append)
    
    return printable_list


def __main__():
    
    # This function governs all other functions and runs
    # only when the program is run directly.
    
    # Read a list of y genes into memory.
    ylist = read_ylist('y_gene_list.txt')
    ylist = list(set(ylist))

    # Read in a dictionary that'll translate Agilent IDs to gene names.
    converter = read_converter('PSG-37 US45103093_251941510006_S01_GE1-v5_95_Feb07_1_4.txt')
    
    # This is the datafile for Zhennan's HSGs. Read that into memory.
    data = read_datafile('2016-08-04_hSG0036_dataMatrix.txt')
    PSG101 = pd.read_csv('2016-03-06_normalized_pt_control_expression_data(1).txt',sep="\t",usecols=['ProbeName','PSG-101 US45103093_251941510019_S01_GE1-v5_95_Feb07_1_2.txt:gProcessedSignal(normalized)'],index_col=0)
    for agilentID in data:
        if agilentID in PSG101.index:
            data[agilentID].append(PSG101.get_value(index=agilentID,col='PSG-101 US45103093_251941510019_S01_GE1-v5_95_Feb07_1_2.txt:gProcessedSignal(normalized)'))
        else:
            data[agilentID].append(0)

    # Get values and ratios for chromosome Y genes from the input data.
    hsg_values = print_file('chrY_gene_values.txt',ylist,converter,data)
    hsg_ratios = print_file('chrY_gene_ratios.txt',ylist,converter,data,ratio=True)
    hsg_dict = {line[0].lower(): line[2:] for line in hsg_ratios}
    
    lower_ylist = [string.lower() for string in np.array(hsg_ratios).transpose()[0][1:]]

    # Here's where the clinical data from our patients comes in.
    clin_data = pd.read_csv('2016-03-17_ClinData_hSG034_parsable.csv')
    patient_rna = pd.read_csv('2016-04-04_rna_features_matrix.csv')
    # Initialize new DataFrames.
    rna_sex_normed = pd.DataFrame(clin_data['SEX'])
    rna_sex_values = pd.DataFrame(clin_data['SEX'])
    # Extract ACTB expression for normalization.
    for column in patient_rna:
        if 'actb_' in column.lower():
            normalization = [pow(value,2) for value in list(patient_rna[column])]
    index = range(28)
    # Iteratively adds a normalized column of RNA
    # expression to that new DataFrame.
    for column in patient_rna:
        flag = 0
        for yvalue in lower_ylist:
            if yvalue in column.split('_')[0].lower():
                flag = 1
        if flag == 1:
            name = patient_rna[column].name
            # Comes in both normed and non!
            rna_sex_normed = rna_sex_normed.join(pd.Series([pow(patient_rna[column][i],2)/normalization[i] for i in range(len(patient_rna[column]))],index=index,name=name))
            rna_sex_values = rna_sex_values.join(pd.Series([pow(value,2) for value in patient_rna[column]],index=index,name=name))

    # We now drop the axis because it was just
    # for the index.
    just_data = rna_sex_normed.drop('SEX',axis=1)
    # Make a new DataFrame!
    new_array = pd.DataFrame(rna_sex_normed['SEX'])

    # Extract the final usable list of chrY genes.
    column_list = []
    for column in just_data:
        column_list.append(column.split('_')[0].lower())
    column_list = list(set(column_list))
    canonical_list = []
    for column in column_list:
        if column in lower_ylist:
            canonical_list.append(column)

    # Iteratively add the chrY genes we care
    # about to our new array, keeping that
    # as the focus of our machine learning.
    names = ['SAMPLE 5','SAMPLE 6','SAMPLE 7','SAMPLE 8','SAMPLE 1','SAMPLE 2','SAMPLE 3','SAMPLE 4','SAMPLE 9','SAMPLE 10','SAMPLE 11','SAMPLE 12','PSG101']
    samples = [names]
    for column1 in canonical_list:
        samples.append(hsg_dict[column1])
        columns_to_check = []
        for column2 in just_data:
            if column2.split('_')[0].lower() == column1:
                columns_to_check.append(just_data[column2].values.tolist())
        columns_to_check = np.array(columns_to_check).transpose()
        final_column = [np.mean(array) for array in columns_to_check]
        new_array = new_array.join(pd.Series(final_column,index=index,name=column1))
    samplesT = np.array(samples).transpose()

    # Split the data into the male and female
    # patients.
    male_array = new_array[new_array.SEX == 'MALE']
    female_array = new_array[new_array.SEX == 'FEMALE']

    # Plot summary statistics of this data.
    length = len(canonical_list)
    i = 0
    for gene in canonical_list:
        i += 1
        plt.clf()
        mzmean = np.mean(male_array[gene])
        female_mean = np.mean(female_array[gene])
        male_std = np.std(male_array[gene])
        female_std = np.std(female_array[gene])
        #plt.axhspan(male_mean-male_std,male_mean+male_std,color='blue',alpha=0.3)
        #plt.axhspan(female_mean-female_std,female_mean+female_std,color='red',alpha=0.3)
        plt.axhspan(male_mean-.5,male_mean+.5,color='blue',alpha=0.5)
        plt.axhspan(female_mean-.5,female_mean+.5,color='red',alpha=0.5)
        plt.scatter(range(13),samples[i])
        plt.xticks(range(13),samples[0],rotation='vertical')
        plt.savefig('gene_plots_ACTB/2017-01-04_'+gene+'_expression_relative_to_sex.pdf',format='pdf',dpi=600)

    # Initialize literally any sklearn model
    # you'd like to (any classifier will work).
    gnb = ske.RandomForestClassifier()
    # Train the model on male/female patients.
    sex_predict = gnb.fit(new_array.drop('SEX',axis=1),new_array['SEX'])

    # Just to get an idea, score yourself.
    # (Might be overfit but very rough idea.)
    print "Score of the predictor is",sex_predict.score(new_array.drop('SEX',axis=1),new_array['SEX'])

    # Now, test the predictor on your samples.
    for sample in samplesT:
        sample_array = np.array([float(value) for value in sample[1:]]).reshape(1,-1)
        print "The prediction for",sample[0],"is",sex_predict.predict(sample_array)


# The condition that the program is run directly.
# This makes unit testing easy.
if __name__ == '__main__':
    __main__()









