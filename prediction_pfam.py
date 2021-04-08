import pandas as pd
import progressbar
import gc
import sys


def pfam_parse(filename):
    rlist = []
    try:
        with open(filename) as file:
            file.seek(0, 0)
            for line in progressbar.progressbar(file):
                if line.startswith("#"):
                    continue
                temp_dict = {}
                x = [y for y in line.split(" ") if y != '']
                if x[9] != "1":
                    continue
                temp_dict["gene_stable_id"] = x[3]
                temp_dict["accession"] = x[1]
                temp_dict["tlen"] = x[2]
                temp_dict["qlen"] = x[5]
                temp_dict["domain"] = x[0]
                temp_dict["hmm_from"] = x[15]
                temp_dict["hmm_to"] = x[16]
                temp_dict["ali_from"] = x[17]
                temp_dict["ali_to"] = x[18]
                temp_dict["env_from"] = x[19]
                temp_dict["env_to"] = x[20]
                rlist.append(temp_dict)
    except Exception as e:
        print(e)
        return pd.DataFrame()
    print(len(rlist))
    tdf = pd.DataFrame(rlist)
    return tdf


def main():
    arg = sys.argv
    fname_1 = arg[-2]
    fname_2 = arg[-1]
    df = pfam_parse(fname_1)
    df.to_hdf("pfam_db_positive.h5", key="pfam_db_positive", mode="w")
    df = ""
    gc.collect()
    df = pfam_parse(fname_2)
    df.to_hdf("pfam_db_negative.h5", key="pfam_db_negative", mode="w")
    print("Pfam Databases Written Successfully :)")


if __name__ == "__main__":
    main()

