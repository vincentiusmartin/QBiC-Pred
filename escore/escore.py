
import time
import os
import gzip
import pandas as pd

chrpath = "/Users/vincentiusmartin/Research/MutationPredictor/QBiC-Pred/chromosomes"

def get_chromosome(cidx,hgver):
    with gzip.open("%s/%s/chr.%s.fa.gz" % (chrpath,hgver,cidx,),'rb') as f:
        next(f)
        chrom = f.read().decode('utf-8').replace('\n','')
    return chrom

def get_wildmut_from_chr(pos,kmer,cidx,hgver):
    with gzip.open("%s/%s/chr.%s.fa.gz" % (chrpath,hgver,cidx,),'rb') as f:
        next(f)
        chrom = f.read().decode('utf-8').replace('\n','')
    seq_wild = chromosome[pos-kmer+1:pos+kmer]  #-5,+6
    seq_mut = chromosome[pos-kmer+1:pos] + row['mutated_to_allele'] + chromosome[pos + 1:pos + kmer]

    return [seq_wild,seq_mut]

def read_intbl(intbl_path,kmer):
    tsv = pd.read_csv(intbl_path,
            sep="\t",
            usecols=['chromosome','chromosome_start','mutation_type','mutated_from_allele','mutated_to_allele'])
    tsv = tsv[tsv['mutation_type'].apply(lambda x: "single base substitution" == x)].drop('mutation_type',1).drop_duplicates() # only take single base mutation
    grouped = tsv.groupby('chromosome',sort=True)
    dataset = {key:item for key,item in grouped}

    result = list()
    for cidx in [str(a) for a in range(1,23)] + ['X','Y']:
        if cidx not in dataset:
            continue
        print("Iterating dataset for chromosome {}...".format(cidx))
        chromosome = get_chromosome(cidx,"hg19")
        for idx,row in dataset[cidx].iterrows():
            pos = row['chromosome_start'] - 1
            if row['mutated_from_allele'] != chromosome[pos]:
                sys.exit("Found mismatch in the mutation: \n{}".format(row))
            seq_wild = chromosome[pos-kmer+1:pos+kmer]  #-5,+6
            seq_mut = chromosome[pos-kmer+1:pos] + row['mutated_to_allele'] + chromosome[pos + 1:pos + kmer]
            # for escore, just use 8?
            result.append([idx,seq_wild,seq_mut]) #rowidx,seq,escore_seq,val,diff,t,pbmname
    return result

'''
Convert integer representation of a sequence to
the sequence string
'''
def itoseq(seqint,kmer):
    nucleotides = {0:'A',1:'C',2:'G',3:'T'}
    binrep = 0
    seq = ""
    while(seqint > 0):
        seq = nucleotides[seqint & 3] + seq
        seqint >>= 2
    while len(seq) < kmer:
        seq = 'A' + seq
    return seq

def seqtoi(seq):
    nucleotides = {'A':0,'C':1,'G':2,'T':3}
    binrep = 0
    for i in range(0,len(seq)):
        binrep <<= 2
        binrep |= nucleotides[seq[i]]
    return binrep

def isbound_escore(seq,etable,kmer=8):
    bsite_cutoff = 0.4
    nbsite_cutoff = 0.3
    nucleotides = {'A':0,'C':1,'G':2,'T':3}
    grapper = (2<<(8*2-1))-1
    binrep = seqtoi(seq[0:kmer])
    elist = [etable[binrep]]
    for i in range(kmer,len(seq)):
        binrep = ((binrep << 2) | seqtoi(seq[i])) & grapper
        elist.append(etable[binrep])
    if max(elist) < nbsite_cutoff:
        return "unbound"
    else:
        isbound = False
        for i in range(0,len(elist)):
            if elist[i] > bsite_cutoff:
                if isbound:
                    return "bound"
                else:
                    isbound = True
            else:
                isbound = False
        return "ambiguous"

def calculate_escore(seq18mer,pbm_name):
    eshort_path = "Mus_musculus_M00423_1.94d_Zoo_01_1289_8mers.txt"
    emap_path = "index_short_to_long.csv"

    # maybe use database?
    with open(eshort_path) as f:
        eshort = [float(line) for line in f]
    with open(emap_path) as f:
        next(f)
        emap = [int(line.split(",")[1])-1 for line in f]

    elong = [eshort[idx] for idx in emap]

    print(seq18mer)
    wild = seq18mer[:-1]
    mut = seq18mer[:8] + seq18mer[-1] + seq18mer[9:-1]

    return isbound_escore(wild,elong),isbound_escore(wild,elong)

'''
    eshort_path = "Mus_musculus_M00423_1.94d_Zoo_01_1289_8mers.txt"
    emap_path = "index_short_to_long.csv"
    mut_data = "JunD_standard_query.txt"
'''
def escore_from_mutdata(eshort_path,emap_path,mut_data):
    #print(calculate_escore("CCCCCATCCAGGGGCCGT",""))

    s = time.time()

    with open(eshort_path) as f:
        eshort = [float(line) for line in f]
    with open(emap_path) as f:
        next(f)
        emap = [int(line.split(",")[1])-1 for line in f]

    elong = [eshort[idx] for idx in emap]

    with open(mut_data) as f:
        mutquery = [line.strip().split("\t") for line in f]

    csvstr = "ref_seq,ref_stats,mut_seq,mut_stats\n"
    for i in range(0,len(mutquery)):
        refseq = mutquery[i][0]
        mutseq = mutquery[i][1]
        csvstr += "%s,%s,%s,%s\n" % (refseq,isbound_escore(refseq,elong),mutseq,isbound_escore(mutseq,elong))

    with open("bounded.csv",'w') as f:
        f.write(csvstr[:-1])

    print("Elapsed time: %.2f" % (time.time()-s))

def escore_from_pbm(seq,pbm_escore):
    emap_path = "escore/index_short_to_long.csv"

    with open(pbm_escore) as f:
        pbm_eshort = [float(line) for line in f]
    with open(emap_path) as f:
        next(f)
        emap = [int(line.split(",")[1])-1 for line in f]

    elong = [pbm_eshort[idx] for idx in emap]

    print(isbound_escore(seq[0],elong),isbound_escore(seq[1],elong))

def update_scorename(escore_dir):
    for file in os.listdir(escore_dir):
        newfilename = file.replace("_8mers","_escore")
        os.rename(escore_dir + "/" + file, escore_dir + "/" + newfilename)


if __name__=="__main__":
    kmer = 8
    #ch = get_chromosome(7,"hg19")
    # ,"/Users/vincentiusmartin/Research/MutationPredictor/QBiC-Pred/QBiC-Pred/mapping_data/pbmtohugo.txt"
    #tbl = read_intbl("../mut-short.tsv",6)
    inlist = [[24, 'CCATCCAGGGG', 'CCATCTAGGGG'], [57, 'CACCTTGGCCT', 'CACCTCGGCCT'], [94, 'GCATGGTGCGG', 'GCATGATGCGG'], [124, 'ACAAAGATACT', 'ACAAAAATACT'], [151, 'GTAAAAGTGAT', 'GTAAACGTGAT'], [191, 'AAAGTGGTATT', 'AAAGTAGTATT'], [246, 'TATGCTGAAAT', 'TATGCGGAAAT'], [262, 'GAACTAGAAAC', 'GAACTGGAAAC'], [313, 'CTACTTTTTGT', 'CTACTGTTTGT'], [314, 'TACCTCATTTA', 'TACCTTATTTA'], [315, 'TCTACCACATT', 'TCTACTACATT'], [334, 'GTGCTGAGATT', 'GTGCTAAGATT'], [406, 'CCACATAATCA', 'CCACAAAATCA'], [426, 'CAAGGCGCGTG', 'CAAGGTGCGTG'], [439, 'GTGTTCCCTGC', 'GTGTTTCCTGC'], [456, 'CCTCCCGATGA', 'CCTCCAGATGA'], [482, 'TTGGTCCACTG', 'TTGGTGCACTG'], [494, 'GAGGCGCAGGA', 'GAGGCACAGGA'], [518, 'GGGGAGGGGGG', 'GGGGACGGGGG'], [608, 'AGAGGAAAAGA', 'AGAGGTAAAGA']]
    escore_from_pbm(['CCATCCAGGGG', 'CCATCTAGGGG'],"escore_rapid/Homo_sapiens|M01221_1.94d|Barrera2016|VSX1_R166Q_R1_escore.txt")
