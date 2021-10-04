from ftplib import FTP
import progressbar
import os
import sys
import requests
import numpy as np

print("Working directory:")
print(os.getcwd())
print("\n")
print("Make download_links directory if necessary")
os.makedirs("download_links", exist_ok=True)

# This wil download all the fasta files for the coding sequences. To change the directory, change the argument in the get_data_file argument.
host = "ftp.ensembl.org"
user = "anonymous"
password = ""

print("Connecting to {}".format(host))
ftp = FTP(host)
ftp.login(user, password)
print("Connected to {}".format(host))
base_link = "ftp://ftp.ensembl.org"
# find sequences of all the cds files
l = ftp.nlst("/pub/release-96/fasta")
lt = []
for x in progressbar.progressbar(l):
    y = ftp.nlst(x + "/cds")
    for z in y:
        if z.endswith(".cds.all.fa.gz"):
            lt.append(z)

print("seq files links obtained.")
# write the file out using numpy
np.savetxt("seq_link.txt", np.array(lt), fmt="%s")
print("seq_link.txt written")

# find all the files with protein sequences
print("Find all fasta files")
l = ftp.nlst("/pub/release-96/fasta")
lt = []
print(type(l))
for x in progressbar.progressbar(l):
    y = ftp.nlst(x + "/pep")
    for z in y:
        if z.endswith(".pep.all.fa.gz"):
            lt.append(z)

print("protein files links obtained.")
# save the protein text file
np.savetxt("protein_seq.txt", np.array(lt), fmt="%s")
print("protein_seq.txt written")

# get link of all the gtf files
l = ftp.nlst("/pub/release-96/gtf")
lt = []
for x in progressbar.progressbar(l):
    y = ftp.nlst(x)
    for z in y:
        if z.endswith(".96.gtf.gz"):
            lt.append(z)

print("gtf files links obtained.")
# save the gtf link files
np.savetxt("gtf_link.txt", np.array(lt), fmt="%s")
print("gtf_link.txt written")

