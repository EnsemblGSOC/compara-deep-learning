import json
import requests
import progressbar
import sys


def update_protein(gene_seq, gene):
    t = 0
    while(t != 2):
        try:
            server = "https://rest.ensembl.org"
            ext = "/sequence/id/" + \
                str(gene) + "?type=protein;multiple_sequences=1"

            r = requests.get(
                server + ext,
                headers={
                    "Content-Type": "application/json"})

            if not r.ok:
                r.raise_for_status()
                sys.exit()
            r = r.json()
            if len(r) == 1:
                r = dict(r[0])
                gene_seq[gene] = str(r["seq"])
                return
            else:
                maxi = 0
                maxlen = 0
                for i in range(len(r)):
                    m = r[i]
                    m = dict(m)
                    if len(m["seq"]) > maxlen:
                        maxi = i
                r = dict(r[maxi])
                gene_seq[gene] = str(r["seq"])
                return
        except BaseException:
            t += 1
            # print("\nError:",e)
            continue
    gene_seq[gene] = ""


def update_rest_protein(data):
    gids = {}
    with open("processed/not_found.json", "r") as file:
        gids = dict(json.load(file))

    gids = list(gids.keys())

    geneseq = {}

    server = "https://rest.ensembl.org"
    ext = "/sequence/id?type=protein"
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"}

    for i in progressbar.progressbar(range(0, len(gids) - 50, 50)):
        ids = dict(ids=list(gids[i:i + 50]))
        while(1):
            try:
                r = requests.post(
                    server + ext,
                    headers=headers,
                    data=str(
                        json.dumps(ids)))
                if not r.ok:
                    r.raise_for_status()
                gs = r.json()
                tgs = {}
                for g in gs:
                    tgs[g["query"]] = g["seq"]
                geneseq.update(tgs)
                break
            except Exception as e:
                print("Error:", e)
                continue

    data.update(geneseq)
    for genes in gids:
        try:
            _ = data[genes]
        except BaseException:
            print(genes)
            update_protein(data, genes)

    print("Gene Sequences Updated Successfully")
    return data


def update(gene_seq, gene):
    t = 0
    while(t != 2):
        try:
            server = "https://rest.ensembl.org"
            ext = "/sequence/id/" + \
                str(gene) + "?type=cds;multiple_sequences=1"

            r = requests.get(
                server + ext,
                headers={
                    "Content-Type": "application/json"})

            if not r.ok:
                r.raise_for_status()
                sys.exit()
            r = r.json()
            if len(r) == 1:
                r = dict(r[0])
                gene_seq[gene] = str(r["seq"])
                return
            else:
                maxi = 0
                maxlen = 0
                for i in range(len(r)):
                    m = r[i]
                    m = dict(m)
                    if len(m["seq"]) > maxlen:
                        maxi = i
                r = dict(r[maxi])
                gene_seq[gene] = str(r["seq"])
                return
        except BaseException:
            t += 1
            # print("\nError:",e)
            continue
    gene_seq[gene] = ""


def update_rest(data, fname):
    gids = {}
    with open("processed/not_found_" + fname + ".json", "r") as file:
        gids = dict(json.load(file))

    gids = list(gids.keys())

    geneseq = {}

    server = "https://rest.ensembl.org"
    ext = "/sequence/id?type=cds"
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"}

    for i in progressbar.progressbar(range(0, len(gids) - 50, 50)):
        ids = dict(ids=list(gids[i:i + 50]))
        while(1):
            try:
                r = requests.post(
                    server + ext,
                    headers=headers,
                    data=str(
                        json.dumps(ids)))
                if not r.ok:
                    r.raise_for_status()
                gs = r.json()
                tgs = {}
                for g in gs:
                    tgs[g["query"]] = g["seq"]
                geneseq.update(tgs)
                break
            except Exception as e:
                print("Error:", e)
                continue

    data.update(geneseq)
    for genes in gids:
        try:
            _ = data[genes]
        except BaseException:
            print(genes)
            update(data, genes)

    print("Gene Sequences Updated Successfully")
    return data
