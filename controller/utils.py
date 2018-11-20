import gzip

def get_chrom(cfile):
    with gzip.open(cfile,'rb') as f:
        next(f)
        chrom = f.read().decode('utf-8').replace('\n','')
    return chrom

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

'''
does not append 1, used for integer indexing
'''
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

"""
return: "is bound wild > is bound mut"
"""
def isbound_escore_18mer(seq18mer,pbm_name):
    eshort_path = "%s/%s_escore.txt" % (app.config['ESCORE_DIR'],pbm_name)
    # TODO: avoid IO, maybe using global var?
    short2long_map = "%s/index_short_to_long.csv" % (app.config['ESCORE_DIR'])

    #  -- this definitely needs to go to a database
    with open(eshort_path) as f:
        eshort = [float(line) for line in f]
    with open(short2long_map) as f:
        next(f)
        emap = [int(line.split(",")[1])-1 for line in f]

    elong = [eshort[idx] for idx in emap]

    wild = seq18mer[:-1]
    mut = seq18mer[:8] + seq18mer[-1] + seq18mer[9:-1]

    return "%s>%s" % (isbound_escore(wild,elong),isbound_escore(mut,elong))

# https://stackoverflow.com/questions/2130016/splitting-a-list-into-n-parts-of-approximately-equal-length
def chunkify(lst,n):
    return [lst[i::n] for i in range(n)]
