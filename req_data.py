import os
import pandas
import sys
import urllib.request as urllib
import pandas as pd


lf=[]

dir_name="data"

def get_data_file(file):
    print(file)
    if not os.path.isfile(file):
        print("The specified file does not exist!!!")
        sys.exit(1)

    with open(file,"r")as f:
        lf=f.read().splitlines()

    if os.path.exists("data"):
        print("Data Directory Already Exists!!!")

    os.mkdir(dir_name)
    for x in lf:
        fname=x.split("/")[-1]
        path=os.path.join(dir_name,fname)
        urllib.urlretrieve(x,path)

    return dir_name
