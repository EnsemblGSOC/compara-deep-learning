from ftplib import FTP
import progressbar
import os
import urllib.request as urllib


def download_data(x, dir_name):
    fname = x.split("/")[-1]
    path = os.path.join(dir_name, fname)
    urllib.urlretrieve(x, path)


def get_data_file(file, dir):
    if not os.path.isfile(file):
        print("The specified file does not exist!!!")
        sys.exit(1)

    with open(file, "r") as f:
        lf = f.read().splitlines()

    if not os.path.exists(dir):
        os.mkdir(dir)
    for x in progressbar.progressbar(lf):
        download_data(x, dir)


os.mkdir("/downloads")
os.mkdir("downloads/gtfs")
os.mkdir("downloads/cds")
os.mkdir("downloads/pep")

print("Downloading Data")
get_data_file("./gtfs/gtf_link.txt", "data")
get_data_file("./cds/seq_link.txt", "geneseq")
get_data_file("./pep/protein_seq.txt", "pro_seq")
print("Download Complete.................")
