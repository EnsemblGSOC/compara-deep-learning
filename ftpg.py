from ftplib import FTP
from req_data import get_data_file


#This wil download all the fasta files for the coding sequences. To change the directory, change the argument in the get_data_file argument.
host ="ftp.ensembl.org"
user = "anonymous"
password = ""

print("Connecting to {}".format(host))
ftp = FTP(host)
ftp.login(user, password)
print("Connected to {}".format(host))
l=ftp.nlst("/pub/release-96/fasta")
lt=[]
for x in l:
    y=ftp.nlst(x+"/cds")
    for z in y:
        if z.endswith(".cds.all.fa.gz"):
            lt.append(z)

base_link="ftp://ftp.ensembl.org"
with open("seq_link.txt","w") as file:
    for x in lt:
        file.write(base_link+x)
        file.write("\n")

l=ftp.nlst("/pub/release-96/gtf")
lt=[]
for x in l:
    y=ftp.nlst(x)
    for z in y:
        if z.endswith(".96.gtf.gz"):
            lt.append(z)

base_link="ftp://ftp.ensembl.org"
with open("gtf_link.txt","w") as file:
    for x in lt:
        file.write(base_link+x)
        file.write("\n")

ch=input("Do you want to download the data?[y/n]")
if ch=='y':
    print("Downloading Data.................")
    get_data_file("gtf_link.txt","data")
    get_data_file("seq_link.txt","geneseq")
    print("Download Complete.................")
