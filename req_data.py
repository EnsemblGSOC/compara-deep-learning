import os
import pandas
import sys
import urllib.request as urllib
import pandas as pd


lf=[]

def download_data(x,dir_name):
    fname=x.split("/")[-1]
    path=os.path.join(dir_name,fname)
    urllib.urlretrieve(x,path)

def get_data_file(file,dir):
    if not os.path.isfile(file):
        print("The specified file does not exist!!!")
        sys.exit(1)

    with open(file,"r")as f:
        lf=f.read().splitlines()

    if not os.path.exists(dir):
        os.mkdir(dir)
    for x in lf:
        download_data(x,dir)
