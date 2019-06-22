import json
import requests
import progressbar
from create_synteny_matrix import update

def update_rest(data):
    gids={}
    with open("processed/not_found.json","r") as file:
        gids=dict(json.load(file))

    gids=list(gids.keys())

    geneseq={}

    server = "https://rest.ensembl.org"
    ext = "/sequence/id?type=cds"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}

    for i in progressbar.progressbar(range(0,len(gids)-50,50)):
        ids=dict(ids=list(gids[i:i+50]))
        while(1):
            try:
                r = requests.post(server+ext, headers=headers, data=str(json.dumps(ids)))
                if not r.ok:
                    r.raise_for_status()
                gs=r.json()
                tgs={}
                for g in gs:
                    tgs[g["query"]]=g["seq"]
                geneseq.update(tgs)
                break
            except Exception as e:
                print("Error:",e)
                continue

    data.update(geneseq)
    for genes in gids:
        try:
            _=data[genes]
        except:
            print(genes)
            update(data,genes)

    with open("processed/gene_sequences.json","w") as file:
        json.dump(data,file)

    print("Gene Sequences Updated Successfully")
    return data
            
