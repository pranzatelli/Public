# Dear Mr. Michael,
# Thank you so much for doing the hard SQL work!
# It really means a lot.
#
# Thank you again,
# Thomas Pranzatelli

import csv
import sys
csv.field_size_limit(sys.maxsize)

import json
import enst_identifier as ei
import copy
import argparse
import os


def get_size(obj, seen=None):
    """Recursively finds size of objects"""
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    return size


def parse_arguments():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--project',type=str)
    parser.add_argument('-g','--memory',type=int)
    
    return parser.parse_args()


def load_json(filepath):
    
    with open(filepath) as file:
        return json.load(file)

def package_as_json(object,number,path,encsr):

    with open(path+'/'+encsr+'-subdictionary'+str(number)+'.json','w') as infile:
        print infile
        json.dump(object,infile)


def dict_to_json(dma,pfp,efp,mem,project,encsr):

    size_before = 0
    count = 0
    dict = {}
    for ENSG in dma:
        ENST_index = dma[ENSG][0]
        enhancer_index = dma[ENSG][1]
        newdict = {}
        for ENST in ENST_index:
            promoter = pfp[ENST]
            ENST_key = promoter[0][-1] + '@' + promoter[0][0]+':'+str(promoter[0][1])+'-'+str(promoter[0][2])
            fp_dict = {}
            for footprint in promoter[1]:
                PWM_dict = {'p-value':footprint[1],'q-value':footprint[2],'score':footprint[3]}
                fp_key = footprint[0][0]+':'+str(footprint[0][1])+'-'+str(footprint[0][2])
                for PWM in footprint[4]:
                    PWM_key = PWM[4]+'@'+PWM[0]+':'+str(PWM[1])+'-'+str(PWM[2])
                    PWM_dict[PWM_key] = PWM
                fp_dict[fp_key] = PWM_dict
            newdict[ENST_key] = fp_dict
        for enhancer_i in enhancer_index:
            enhancer = efp[enhancer_i]
            enhancer_key = enhancer[0][0]+':'+str(enhancer[0][1])+'-'+str(enhancer[0][2])
            fp_dict = {}
            for footprint in enhancer[1]:
                PWM_dict = {'p-value':footprint[1],'q-value':footprint[2],'score':footprint[3]}
                fp_key = footprint[0][0]+':'+str(footprint[0][1])+'-'+str(footprint[0][2])
                for PWM in footprint[4]:
                    PWM_key = PWM[4]+'@'+PWM[0]+':'+str(PWM[1])+'-'+str(PWM[2])
                    PWM_dict[PWM_key] = PWM
                fp_dict[fp_key] = PWM_dict
            newdict[enhancer_key] = fp_dict
        size_after = get_size(newdict)
        if size_after + size_before >= mem:
            package_as_json(dict,count,project,encsr)
            dict = {}
            size_before = 0
            count += 1
        size_before = size_before+size_after
        dict[ENSG] = newdict
    package_as_json(dict,count,project,encsr)


def __main__():
    
    args = parse_arguments()
    project = args.project
    memory = args.memory * 750000000
    encsr = project.split('/')[-1]
    
    e_fp = load_json(project+'/'+encsr+'-enhancer_footprints.json')
    p_fp = load_json(project+'/'+encsr+'-promoter_footprints.json')
    ei.append_enst(p_fp)
    d_ma = load_json(project+'/'+encsr+'-meta_associations.json')

    dict_to_json(d_ma,p_fp,e_fp,memory,project,encsr)



if __name__ == '__main__':
    __main__()
