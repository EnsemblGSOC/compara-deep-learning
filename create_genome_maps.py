import pickle
from get_data import get_data_genome

dir_g = "data"
cmap, cimap, ld, ldg, a, d = get_data_genome(dir_g)

data = dict(cmap=cmap, cimap=cimap, ld=ld, ldg=ldg, a=a, d=d)

with open("genome_maps", "wb") as file:
    pickle.dump(data, file)

print("Genome Maps Created Successfully.")
