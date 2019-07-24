from read_data import read_data_genome
from process_data import list_dict_genomes, create_chromosome_maps


def get_data_genome(dir):
    a = []
    d = {}
    ld = []
    ldg = []
    a, d = read_data_genome(dir, a, d)
    assert(len(a) == len(d))
    print("Creating Maps:")
    ld, ldg = list_dict_genomes(a, d)
    cmap, cimap = create_chromosome_maps(a, d)
    assert(len(ld) == len(ldg))
    for i in range(len(ld)):
        assert(len(ld[i]) == len(ldg[i]))
    return cmap, cimap, ld, ldg, a, d
