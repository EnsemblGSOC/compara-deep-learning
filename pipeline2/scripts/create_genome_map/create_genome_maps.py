import pickle
from get_data import get_data_genome
import sys

# read in the input directory
dir_g = sys.argv[1]
out_path = sys.argv[2]

cmap, cimap, ld, ldg, a, d = get_data_genome(dir_g)

data = dict(cmap=cmap, cimap=cimap, ld=ld, ldg=ldg, a=a, d=d)

with open(out_path, "wb") as file:
    pickle.dump(data, file)

print("Genome Maps Created Successfully.")
