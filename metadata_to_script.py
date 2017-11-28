import csv
import sys
csv.field_size_limit(sys.maxsize)

import os


def read_metadata_to_memory():
    
    # This is a function that reads a .csv file
    # of the metadata and relays that
    # to a list object.
    
    # Here we open the .csv file.
    metadata_csv = open('/data/ChioriniCompCor/metamachine/taskoffice/metadata.csv','rU')
    
    # Here we initialize the returned list.
    return_list = []
    
    # Here we iterate through the .csv file,
    # pulling each comma-delimited line as a
    # list into our list.
    for line in csv.reader(metadata_csv, delimiter='\n'):
        for comma in csv.reader(line, delimiter=','):
            return_list.append(comma)

    # Closing files is hygienic.
    metadata_csv.close()
    
    # Remove the header and return the list.
    return return_list[1:]


def check_script_number():
    
    # The purpose of this function is to find the
    # highest number given to a script thus far,
    # and name scripts with numbers starting at
    # that number.
    
    # Here we read in the names all the scripts
    # that have and haven't been used in the
    # taskoffice.
    scripts = os.listdir('/data/ChioriniCompCor/metamachine/taskoffice/scripts')
    donescripts = os.listdir('/data/ChioriniCompCor/metamachine/taskoffice/donescripts')
    
    # It's important to convert these numbers to
    # integers so that you can find the maximum.
    numbers = [int(number.split('.')[0]) for number in scripts + donescripts]
    maxnumber = max(numbers)

    return maxnumber


def check_data_folder():
    
    # The purpose of this function is to make a
    # list of the ENCODE directories so we don't
    # write scripts to produce something that
    # already exists.
    
    # Here we read in the names all the directories
    # in the ENCODE data directory.
    ENCODE = os.listdir('/data/ChioriniCompCor/metamachine/ENCODE')
    
    return ENCODE


def dict_from_metadata(metadata):

    # The purpose of this function is to take the
    # large, messy metadata structure and turn it
    # into a more parseable dict for later ease.
    
    # Here we initialize the returned list.
    return_dict = {}
    
    # Here we iterate over the function.
    # line[6] is the ENCSR project, and the key
    # used in our dict. At that key, the dict
    # returns a list of three lists: the ENCFFs
    # for PE1, PE2 and SE reads.
    for line in metadata:
        if line[6] not in return_dict:
            return_dict[line[6]] = [[],[],[]]
        if 'P' in line[10]:
            if line[11] == '1':
                return_dict[line[6]][0].append(line[8])
            if line[11] == '2':
                return_dict[line[6]][1].append(line[8])
        if 'S' in line[10]:
            return_dict[line[6]][2].append(line[8])

    return return_dict


def filter_dict(dict,list):
    
    # The purpose of this function is to take a
    # dictionary and filter out keys based on
    # a given list
    
    # Here we iterate over the function.
    # line[6] is the ENCSR project, and the key
    # used in our dict. At that key, the dict
    # returns a list of three lists: the ENCFFs
    # for PE1, PE2 and SE reads.
    for line in list:
        try: del dict[line]
        except KeyError: pass

    return dict


def write_scripts(dict,number):
    
    # This is a function that turns the final dict
    # into a series of scripts.
    
    # Here we iterate through the dict.
    i = number
    for key in dict:
        i += 1
        # It's important to open and close files
        # to promote hygiene.
        scriptfile = open('/data/ChioriniCompCor/metamachine/taskoffice/scripts/'+str(i)+'.sh','w')
        
        script_to_write = '#!/bin/sh\ncd /data/ChioriniCompCor/metamachine/ENCODE/;mkdir '+key+';cd '+key+';'
        if len(dict[key][0]) > 0:
            script_to_write += 'mkdir P1 P2; cd P1; wget '
            for ENCFF in dict[key][0]:
                script_to_write += 'https://www.encodeproject.org/files/'+ENCFF+'/@@download/'+ENCFF+'.fastq.gz' + ' '
            script_to_write += '; cat * > '+key+'.fastq.gz; rm ENCFF*.fastq.gz ; cd ../P2; wget '
            for ENCFF in dict[key][1]:
                script_to_write += 'https://www.encodeproject.org/files/'+ENCFF+'/@@download/'+ENCFF+'.fastq.gz' + ' '
            script_to_write += '; cat * > '+key+'.fastq.gz; rm ENCFF*.fastq.gz ;cd ..; wget '
        else:
            script_to_write += 'wget '
        for ENCFF in dict[key][2]:
            script_to_write += 'https://www.encodeproject.org/files/'+ENCFF+'/@@download/'+ENCFF+'.fastq.gz' + ' '
        script_to_write += '; cat * > '+key+'.fastq.gz; rm ENCFF*.fastq.gz '
        scriptfile.write(script_to_write)

        scriptfile.close()


def __main__():
    
    # This function governs all other functions and runs
    # only when the program is run directly.
    
    # Here we read the input .tsv file and turn it into a
    # parseable list.
    metadata = read_metadata_to_memory()
    
    # Here we check the highest script number so we can
    # name our scripts later.
    number = check_script_number()

    # Here we check the ENCODE directory to see what
    # hasn't already been downloaded.
    ENCODElist = check_data_folder()

    # Here we iterate over the metadata, producing
    # a list of all the ENCSR projects and the names
    # of their constituent files.
    scriptable_dict = dict_from_metadata(metadata)

    # Here we filter the to-script list and remove
    # elements for which there is already a data file.
    filter_dict(scriptable_dict,ENCODElist)

    # The final step is to write the scripts!
    write_scripts(scriptable_dict,number)


# The condition that the program is run directly.
# This makes unit testing easy.
if __name__ == '__main__':
    __main__()









