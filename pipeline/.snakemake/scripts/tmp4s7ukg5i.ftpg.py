
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/home/aidan/anaconda3/envs/compara/lib/python3.9/site-packages', '/mnt/c/Users/aidan/Google Drive/compara-deep-learning/pipeline/scripts']); import pickle; snakemake = pickle.loads(b'\x80\x04\x95\xe4\x03\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94\x8c\x0fscripts/ftpg.py\x94a}\x94(\x8c\x06_names\x94}\x94\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x10\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h\x16)}\x94\x8c\x05_name\x94h\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94(\x8c\x0cgtf_link.txt\x94\x8c\x0cseq_link.txt\x94\x8c\x0fprotein_seq.txt\x94e}\x94(h\x0c}\x94h\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01e}\x94(h\x0c}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94uh\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bhWK\x01hYK\x01ub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x0e]\x94(h\x10h\x11eh\x10h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x10sNt\x94bh\x11h\x14h\x16\x85\x94R\x94(h\x16)}\x94h\x1ah\x11sNt\x94bub\x8c\x06config\x94}\x94\x8c\x04rule\x94\x8c\x12get_download_links\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8cF/mnt/c/Users/aidan/Google Drive/compara-deep-learning/pipeline/scripts\x94ub.'); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/mnt/c/Users/aidan/Google Drive/compara-deep-learning/pipeline/scripts/ftpg.py';
######## snakemake preamble end #########
from ftplib import FTP
import progressbar
import os
import sys
import requests
import numpy as np


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

