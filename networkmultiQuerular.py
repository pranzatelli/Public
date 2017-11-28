import argparse
import json
import networkx as nx

G = nx.DiGraph()

def parse_arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--project',type=str)
    parser.add_argument('-g','--genelist',type=str)
    parser.add_argument('-d','--depth',type=int)
    parser.add_argument('-n','--name',type=str)

    return parser.parse_args()

def load_json(filepath):

    with open(filepath) as file:
        return json.load(file)

def read_ensg_name(genome):

    openfile = open('/data/ChioriniCompCor/metamachine/genomes/'+genome+'/name_to_ensg.txt')
    readfile = openfile.read()

    outputdict = {}
    for line in readfile.split('\n')[:-1]:
        tabs = line.split('\t')
        outputdict[tabs[1]] = tabs[0]

    return outputdict

def recursive_pwm_search(dict,plist,elist,gene,depth,ensg_name):
    
    if depth == 0:
        return None
    depth = depth - 1
    
    new_genes = []
    if gene in dict:
        p_indices = dict[gene][0]
        e_indices = dict[gene][1]
        for index in p_indices:
            for footprint in plist[index][1]:
                new_genes.append(footprint[1][4])
        for index in e_indices:
            for footprint in elist[index][1]:
                new_genes.append(footprint[1][4])
    
    new_genes = list(set(new_genes))

    if gene in ensg_name:
        gene_proper_name = ensg_name[gene]
    else:
        gene_proper_name = gene

    for new_gene in new_genes:
        recursive_pwm_search(dict,plist,elist,new_gene,depth,ensg_name)
        if new_gene in ensg_name:
            new_gene_proper_name = ensg_name[new_gene]
        else:
            new_gene_proper_name = new_gene
        G.add_edge(new_gene_proper_name,gene_proper_name)


def __main__():

    args = parse_arguments()
    project = args.project
    genes = args.genelist.split(',')
    depth = args.depth
    name = args.name
    
    genome = project.split('_')[0]
    ensg_name = read_ensg_name(genome)

    enhancer_fp = load_json('/data/ChioriniCompCor/metamachine/output-v2/'+project+'-enhancer_footprints.json')
    promoter_fp = load_json('/data/ChioriniCompCor/metamachine/output-v2/'+project+'-promoter_footprints.json')
    index_dict = load_json('/data/ChioriniCompCor/metamachine/output-v2/'+project+'-meta_associations.json')
    for gene in genes:
        recursive_pwm_search(index_dict,promoter_fp,enhancer_fp,gene,depth,ensg_name)

    nx.write_graphml(G,name+'.graphml')

if __name__ == '__main__':
    __main__()
