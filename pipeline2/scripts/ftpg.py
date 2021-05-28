from ftplib import FTP
import progressbar
import os
import sys
import requests
import pandas as pd

print("Working directory:")
print(os.getcwd())
print("\n")
print("Make download_links directory if necessary")
os.makedirs("download_links", exist_ok=True)

# This wil download all the fasta files for the coding sequences. To change the directory, change the argument in the get_data_file argument.
host = "ftp.ensembl.org"
user = "anonymous"
password = ""
# base link for all download paths
base_link = "ftp://ftp.ensembl.org"


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
# save the download links to text file with base_link preprended
ftp_paths = base_link + pd.Series(lt)
ftp_paths.to_csv("seq_link.txt", index=False, header=False)
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
ftp_paths = base_link + pd.Series(lt)
ftp_paths.to_csv("protein_seq.txt", index=False, header=False)
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
ftp_paths = base_link + pd.Series(lt)
ftp_paths.to_csv("gtf_link.txt", index=False, header=False)
print("gtf_link.txt written")

