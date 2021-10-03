import sys
import argparse
import json
import random
import gzip

def load_species_set(file):
    '''
    :param file:
    :return:
    '''
    species_set = []
    with open(file) as file_handle:
        for line in file_handle:
            tab = line.replace("\n","").split()
            species_set.append(tab)
    return species_set

def get_species_from_gene_name(gene, species_set):
    '''
    :param gene:
    :param species_set:
    :return:
    '''

    for sp in species_set:
        if sp[1] in gene:
            return sp[0]
    return None

def belong_to_ss(species, species_set):
    '''
    :param gene:
    :param species_set:
    :return:
    '''
    for sp in species_set:
        if sp[0] == species:
            return True
    return False

def load_families(file, species_set):
    '''
    :param file:
    :param species_set:
    :return:
    '''
    families = {}
    with open(file) as file_handle:
        for line in file_handle:
            tab = line.replace("\n", "").split()
            if len (tab) < 2:
                continue
            if not belong_to_ss(tab[-1], species_set):
                continue
            if not tab[0] in families:
                families[tab[0]] = []
            families[tab[0]].append([tab[1], tab[2]])
    return families

def load_orthologues(orthologues, species_set):
    '''
    :param orthologues:
    :param species_set:
    :return:
    '''
    ortho= []
    print(species_set)
    with gzip.open(orthologues, "rt") as file_handle:
        for line in file_handle:
            # print(type(line))
            tab = line.strip().split("\t")
            # print(tab[2], tab[7])
            if belong_to_ss(tab[2], species_set) and belong_to_ss(tab[7], species_set) and tab[2] != tab[7]:
                # print("WOOH")
                new_tab = tab[:3] + tab[5:8] + tab[4:5]
                # print(new_tab)
                ortho.append(new_tab)
    return ortho

def _get_paralogue(family, ref_sp, number, ss):

    i = 0
    ref_genes = []
    paralogues = []
    while i < len(family):
        gene = family[i]
        if get_species_from_gene_name(gene[0], ss) == ref_sp:
            ref_genes.append(gene)
            family = family[:i] + family[i+1:]
            continue
        i = i + 1
    ref_gene = None
    if len(ref_genes) == 1:
        ref_gene = ref_genes[0]
    elif len(ref_genes) > 1:
        indice = random.randint(0, len(ref_genes) - 1)
        ref_gene = ref_genes[indice]
    else:
        return paralogues
    para = family
    if number < len(family):
        para = random.sample(family, number)
    for par in para:
        tab = ref_gene + [ref_sp] + par + [get_species_from_gene_name(par[0],ss)] + ["paralog"]
        paralogues.append(tab)
    return paralogues



def select_paralogues(families, orthologues, ref_species, ss):
    lst_paralogues = []
    for fam in families.keys():
        genes = families[fam]
        paralogues = _get_paralogue(genes, ref_species, 15, ss)
        for para in paralogues:
            if not para[0] in orthologues or not para[3] in orthologues[para[0]]:
                lst_paralogues.append(para)
    return lst_paralogues

def _get_dico_orthologues(orthologues):
    dico ={}
    for orth in orthologues:
        if not orth[0] in dico:
            dico[orth[0]] = []
        dico[orth[0]].append(orth[3])
    return dico

def _get_gene_from_ref_species(genes, ref_species, ss):
    for gene in genes:
        if get_species_from_gene_name(gene[0],ss) == ref_species:
            return gene
    return None

def _get_gene_from_non_ref_species(gene_non_hom, ref_species, nb, ss):
    genes = []
    for gene in gene_non_hom:
        if get_species_from_gene_name(gene[0], ss) != ref_species:
            if len(genes) < nb:
                genes.append(gene)
            else:
                break
    return genes

def select_negatives(families, ref_species, ss):
    i = 0
    family_names = list(families.keys())
    result_negatives = []
    while i < len(family_names):
        genes = families[family_names[i]]
        gene_ref = _get_gene_from_ref_species(genes, ref_species, ss)
        if gene_ref is None:
            i = i + 1
            continue
        ind_family = 0
        if i < len(family_names)/2:
            ind_family = random.randint(i+1, len(family_names) - 1)
        else:
            ind_family = random.randint(0, i - 1)
        gene_non_hom = families[family_names[ind_family]]
        non_homologues = _get_gene_from_non_ref_species(gene_non_hom, ref_species, 10, ss)
        for neg in non_homologues:
            tab = gene_ref + [ref_species] + neg + [get_species_from_gene_name(neg[0], ss)] + ["non_homologue"]
            result_negatives.append(tab)
        i = i + 1
    return result_negatives

def save_training_data(data_set, out_file):
    print(data_set)
    print(len(data_set))
    with open(out_file, "w") as file_handler:
        for i, data in enumerate(data_set):
            if data == None:
                print(i)
            file_handler.write("\t".join(data)+"\n")



def main(args):
    '''
    :param args:
    :return:
    '''

    nb_elements = int(args.size/3)

    #load the species set
    ss = load_species_set(args.ss)

    # load the families
    fam = load_families(args.fam, ss)

    # load orthologues
    orth = load_orthologues(args.ort, ss)
    print(type(orth))
    print(orth)
    print(orth[0][2])
    # get random orthiologues
    # print(orth)
    random_ortho = random.sample(orth, nb_elements)

    #get random paralogs
    dico_ortho = _get_dico_orthologues(orth)
    paralogues = select_paralogues(fam, dico_ortho, orth[0][2], ss)
    random_para = random.sample(paralogues, nb_elements)

    # get non homologue
    non_homologues = select_negatives(fam, orth[0][2], ss)
    random_non_honmologue = random.sample(non_homologues, nb_elements)

    print(len(random_ortho))
    #print(random_ortho)
    print(len(random_para))
    #print(random_para)
    print(len(random_non_honmologue))
    #print(random_non_honmologue)

    training_set = random_ortho + random_para + random_non_honmologue
    save_training_data(training_set, args.out)


parser = argparse.ArgumentParser()
parser.add_argument('--fam', help='family file' , default="")
parser.add_argument('--ss', help='species set file', default="")
parser.add_argument('--ort', help='file with the orthologues')
parser.add_argument('--out',  help='dataset outfile', default="./train_dataset.txt")
parser.add_argument('--size', type=int, help='size of the dataset', default=1000)
args = parser.parse_args()

main(args)

# python create_dataset.py --fam ../hmm_table.tsv --ss ../config/Valid_species_new_way.txt --ort /nfs/production/flicek/ensembl/production/ensemblftp/release-104/tsv/ensembl-compara/homologies/homo_sapiens/Compara.104.protein_default.homologies.tsv.gz