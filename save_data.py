import os
import pickle
import json
import sys

def write_dict_json(name,dir,d):
    if not os.path.exists(dir):
        os.mkdir(dir)

    path=os.path.join(dir,name+".json")
    with open(path,'w') as file:
        json.dump(d,file)

def write_file_multiple_json(name,dir,l):
        if not os.path.exists(dir):
            os.mkdir(dir)

        path=os.path.join(dir,name)
        with open(path,"w") as file:
            for x in l:
                out=json.dumps(x)
                file.write(out)
                file.write('\n')

def save_data_json(name,data):
    with open(name+".json","w")as file:
        json.dump(data,file)
