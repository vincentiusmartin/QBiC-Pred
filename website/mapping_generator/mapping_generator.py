import os
import pandas as pd
import shutil

def dbd2hugo2pbm(gene_dict,tfdb_df):
    dbd_mapping = {}

    for gene in gene_dict.keys():
        dbds = tfdb_df[tfdb_df["HGNC symbol"]==gene]["DBD"].unique().tolist()
        for dbd in dbds: # a TF shouldn't have >1 DBD, but just in case
            if dbd in dbd_mapping:
                dbd_mapping[dbd].append("%s:%s"%(gene,",".join(gene_dict[gene])))
            else:
                dbd_mapping[dbd] = ["%s:%s"%(gene,",".join(gene_dict[gene]))]
    with open("mapping_data/hugotopbm.txt",'w') as f:
        dbds = sorted(dbd_mapping.keys(), key=lambda s: s.lower())
        towrite = ""
        for dbd in dbds:
            liststr = ";".join(dbd_mapping[dbd])
            towrite += ("%s->%s\n" % (dbd,liststr))
        f.write(towrite.strip())

def pbm2hugo(gene_dict):
    pbm_dict = {}
    for gene in gene_dict.keys():
        for tf in gene_dict[gene]:
            if tf in pbm_dict and gene not in pbm_dict[tf]:
                pbm_dict[tf].append(gene)
            else:
                pbm_dict[tf] = [gene]
    pbms = sorted(pbm_dict.keys(), key=lambda s: s.lower())
    with open("mapping_data/tflist.txt",'w') as f:
        f.write("\n".join(gene_dict.keys()))
    with open('mapping_data/pbmlist.txt','w') as f:
        f.write("\n".join(pbms))
    with open("mapping_data/pbmtohugo.txt",'w') as f:
        towrite = ""
        for pbm in pbms:
            towrite += "%s:%s\n"%(pbm,",".join(pbm_dict[pbm]))
        f.write(towrite.strip())

def generate_mapping_webserv(infile,tfdb):
    gene_df = pd.read_csv(infile)
    gene_dict = dict(zip(gene_df['gene'], gene_df['upbm']))
    gene_dict = {gene:[os.path.splitext(tf)[0].strip() for tf in gene_dict[gene].split(";")] for gene in gene_dict.keys()}
    tfdb_df = pd.read_csv(tfdb)

    # Create a directory to write the mapping results
    if os.path.exists("mapping_data") and os.path.isdir("mapping_data"):
        shutil.rmtree("mapping_data")
    os.mkdir("mapping_data")

    dbd2hugo2pbm(gene_dict,tfdb_df)
    pbm2hugo(gene_dict)

if __name__ == "__main__":
    gene2upbm = "../app/static/files/mapping_gene_to_upbm_final.csv"
    tfdb = "resources_from_humantfs/DatabaseExtract_v_1.01.csv"
    generate_mapping_webserv(gene2upbm,tfdb)
