


class GTFAdaptor:

    def __init__(self, file=""):
        self.chr = {}
        self.file = file
        self.genes = {}

    def load_genes(self, btype=""):
        with open(self.file) as file_handler:
            i = 0
            for line in file_handler:
                if line[0] == "#":
                    continue
                columns = line.split("\t")
                if columns[2] != "gene":
                    continue
                fields = columns[-1].split(";")
                biotype = ""
                gene_id = ""
                gene_name = ""
                gene_source = ""
                for field in fields:
                    if "gene_biotype" in field:
                        biotype = field.split()[-1].replace('"','')
                    elif "gene_id" in field:
                        gene_id = field.split()[-1].replace('"', '')
                    elif "gene_name"  in field:
                        gene_name = field.split()[-1].replace('"', '')
                    elif "gene_source" in field:
                        gene_source = field.split()[-1].replace('"', '')
                if btype != "" and biotype != btype:
                    continue
                gene = {}
                gene["chr"] = columns[0]
                gene["start"] = int(columns[3])
                gene["end"] = int(columns[4])
                gene["strand"] = 1
                if columns[6] == "-":
                    gene["strand"] = -1
                gene["stable_id"] = gene_id
                gene["gene_name"] = gene_name
                gene["biotype"] = biotype
                gene["source"] = gene_source

                if not gene["chr"] in self.chr:
                    self.chr[gene["chr"]]=[]
                    i = 0
                self.chr[gene["chr"]].append(gene)
                self.genes[gene["stable_id"]] = [gene, gene["chr"], i]
                i=i+1


    def getNeighbourGenes(self, gene_id, nb_nghbr, sens="both"):

        if not gene_id  in self.genes:
            raise Exception("Error gene stable id not in the gene data set")
        up_genes = [None] * nb_nghbr
        down_genes = [None] * nb_nghbr
        gene = self.genes[gene_id][0]
        chr = self.genes[gene_id][1]
        index = self.genes[gene_id][2]
        genes = self.chr[chr]

        if sens == "both" or sens == "upstream":
            if index <= nb_nghbr:
                up_genes = genes[:index]
                nones = [None] * (nb_nghbr - len(up_genes))
                up_genes = nones + up_genes
            else:
                up_genes = genes[ index - nb_nghbr : index ]

        if sens == "both" or sens == "downstream":
            if index >= (len(genes) - nb_nghbr):
                down_genes = genes[index+1:]
                nones = [None] * (nb_nghbr - len(down_genes))
                down_genes = down_genes + nones
            else:
                down_genes = genes[index+1:index+nb_nghbr+1]

        return up_genes, down_genes







