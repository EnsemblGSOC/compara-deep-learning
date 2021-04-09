from ftplib import FTP
import progressbar
import os
import sys
import urllib.request as urllib


def download_data(x, dir_name):
    fname = x.split("/")[-1]
    path = os.path.join(dir_name, fname)
    urllib.urlretrieve(x, path)


def get_data_file(file, dir):
    if not os.path.isfile(file):
        print("The specified file does not exist!!!")
        sys.exit(1)

    with open(file, "r")as f:
        lf = f.read().splitlines()

    if not os.path.exists(dir):
        os.mkdir(dir)
    for x in progressbar.progressbar(lf):
        download_data(x, dir)


# This wil download all the fasta files for the coding sequences.
# To change the directory, change the argument in the get_data_file argument.
host = "ftp.ensembl.org"
user = "anonymous"
password = ""

print("Connecting to {}".format(host))
ftp = FTP(host)
ftp.login(user, password)
print("Connected to {}".format(host))
base_link = "ftp://ftp.ensembl.org"
# find sequences of all the cds files
list_of_files = ftp.nlst("/pub/release-96/fasta")
lt = []
for x in list_of_files:
    y = ftp.nlst(x+"/cds")
    for z in y:
        if z.endswith(".cds.all.fa.gz"):
            lt.append(z)

with open("seq_link.txt", "w") as file:
    for x in lt:
        file.write(base_link+x)
        file.write("\n")

# find all the files with protein sequences
list_of_files = ftp.nlst("/pub/release-96/fasta")
lt = []
for x in progressbar.progressbar(list_of_files):
    y = ftp.nlst(x+"/pep")
    for z in y:
        if z.endswith(".pep.all.fa.gz"):
            lt.append(z)

with open("protein_seq.txt", "w") as file:
    for x in lt:
        file.write(base_link+x)
        file.write("\n")
# get link of all the gtf files
list_of_files = ftp.nlst("/pub/release-96/gtf")
lt = []
for x in list_of_files:
    y = ftp.nlst(x)
    for z in y:
        if z.endswith(".96.gtf.gz"):
            lt.append(z)

with open("gtf_link.txt", "w") as file:
    for x in lt:
        file.write(base_link+x)
        file.write("\n")

print("Downloading Data")
get_data_file("gtf_link.txt", "data")
get_data_file("seq_link.txt", "geneseq")
get_data_file("protein_seq.txt", "pro_seq")
print("Download Complete.................")
