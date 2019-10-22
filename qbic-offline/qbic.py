
import pandas as pd
import time
import os
import multiprocessing as mp
import numpy as np
import scipy.stats
import argparse

from timeit import default_timer as timer

import utils
import config

def inittbl(filename, chromosome_version, kmer = 6, filetype=""):
    start = time.time()
    file_extension = os.path.splitext(filename)[1]

    # validate inputs
    if file_extension != ".txt" and file_extension != ".vcf" and file_extension != ".tsv" and file_extension != ".csv":
        raise Exception('File extension should be either txt, vcf, tsv, or csv')

    chrs = [name for name in os.listdir(config.CHRDIR) if os.path.isdir(config.CHRDIR + "/" + name)]
    if chromosome_version not in chrs:
        raise Exception('Chromosome version is not available, please use one of the following: %s' % ", ".join(chrs))

    result = []
    error = ""

    type = ""
    if not filetype:
        if file_extension == ".txt":
            type = "customseq"
        elif file_extension == ".vcf":
            type = "vcf"
        else:
            type = "icgc"
    else:
        type = str(filetype)

    # TODO: if fast enough, we can also put error checking in here
    if type == "customseq":
        with open(filename) as f:
            idx = 0
            for line in f:
                line = line.strip() # removing any whitespace
                if not line.strip(): # empty line
                    continue
                elif "\t" in line:
                    line = line.split("\t")
                else:
                    line = line.split()
                idx += 1
                # Input checking
                if line[1] == "*":
                    print("Found * (deletion) in the alt field. Currrently it is not handled by QBiC, this line will be skipped")
                    continue
                elif line[1].upper() not in ['A','C','G','T']:
                    print("%s is not a valid nucleotide (A,C,G,T)")
                    continue
                # line[1] is the base mid nucleotide mutated to
                escore_seq = line[0] + line[1]
                mid_seq = escore_seq[len(escore_seq)//2-6:len(escore_seq)//2+5] + line[1] # the 12mer seq
                result.append([idx,mid_seq,escore_seq,utils.seqtoi(mid_seq),0,0,"None"])
    else:
        if type == "vcf":
            df = pd.read_csv(filename,sep="\t",header=None).drop(2,1)
            df = df.rename(columns={0:"chromosome",1:"pos",3:"mutated_from",4:"mutated_to"})
        elif type == "icgc" or type == "mut":
            if file_extension == ".tsv":
                separator = "\t"
            else: # must be csv since we checked it, TODO: can also return error here
                separator = ","
            if type == "icgc":
                df = pd.read_csv(filename,
                            sep=separator,
                            usecols=['chromosome','chromosome_start','mutation_type','mutated_from_allele','mutated_to_allele'])
                df = df[df['mutation_type'].apply(lambda x: "single base substitution" == x)].drop('mutation_type',1).drop_duplicates() # only take single base mutation
                df = df.rename(columns={"chromosome_start":"pos","mutated_from_allele":"mutated_from","mutated_to_allele":"mutated_to"})
            else: # type == "mut", usually have chr, start, end, ref, alt
                df = pd.read_csv(filename, sep=separator)
                df.columns = ["chromosome",'pos','chromosome_end',"mutated_from","mutated_to"]
                df = df.drop('chromosome_end',1)
        else:
            raise Exception("File type is not recognized")

        # remove "chr" from all chromosomes
        df['chromosome'] = df['chromosome'].map(lambda x:x.replace("chr",""))

        grouped = df.groupby('chromosome',sort=True)
        dataset = {key:item for key,item in grouped}

        for cidx in [str(a) for a in range(1,23)] + ['X','Y']:
            if cidx not in dataset:
                continue
            print("Iterating dataset for chromosome {}...".format(cidx))
            chromosome = utils.get_chrom(config.CHRDIR + "/" + chromosome_version + "/chr." + str(cidx) + '.fa.gz')
            for idx,row in dataset[cidx].iterrows():
                pos = row['pos'] - 1
                if row['mutated_from'] != chromosome[pos]:
                    cver = cpath.split("/")[-1]
                    error = "For the input mutation %s>%s at position %s in chromosome %s, the mutated_from nucleotide (%s) does not match the nucleotide in the %s reference genome (%s). Please check the input data and verify that the correct version of the reference human genome was used." % (row['mutated_from'], row['mutated_to'], row['pos'], row['chromosome'], row['mutated_from'], cver, chromosome[pos])
                    raise Exception(error)
                seq = chromosome[pos-kmer+1:pos+kmer] + row['mutated_to'] #-5,+6
                # for escore, just use 8?
                escore_seq = chromosome[pos-9+1:pos+9] + row['mutated_to']
                result.append([idx,seq,escore_seq,utils.seqtoi(seq),0,0,"None"]) #rowidx,seq,escore_seq,val,diff,t,pbmname

    result = sorted(result,key=lambda result:result[0])
    # example row in result: [73, 'CCAACCAACCCA', 'ATTCCAACCAACCCCCTA', 5263444, 0, 0, 'None']
    print("Time to preprocess: {:.2f}secs".format(time.time()-start))
    return result

def predict(predlist, dataset, ready_count,
            filteropt="p-value", filterval=0.001, spec_ecutoff=0.4, nonspec_ecutoff=0.35):
    """
    for the container list, key is a tuple of: (rowidx,sequence,seqidx)
    and each element in value is a list if: [diff,z-score,pbmname]
    If spec_ecutoff and nonspec_ecutoff == -1 then no escore calculation is done

    filteropt = p-value, t-value, none
                if none then just calculate escore

    return:
     filteropt=1: diff,z_score,tfname
     filteropt=2: diff,p_val,escore,tfname
    """
    buggedtf = 0
    #[96, 'TCATGGTGGGTT', GCTTCATGGTGGGTGGAT, 13872815, 0, 0, '-'] -- 37, 'GCCCAGAAAGGA', 9773096
    if filteropt == "p-value": #t-value
        # leave this empty as for p-value, we don't have to compare and the size is dynamic
        container = {tuple(row[:4]):[] for row in dataset}
    else:
        container = {tuple(row[:4]):[[0,0,1,"None","None"]] for row in dataset} # rowidx,12mer,18mer,seqidx : [diff,z,p,bind,pbmname]

    escore_total_time = 0
    # iterate for each transcription factor
    for i in range(0,len(predlist)):
        start = time.time()
        pbmname = '.'.join(map(str,predlist[i].split(".")[1:-1]))
        print("Processing " + pbmname)
        with open(predlist[i], 'r') as f:
            tflist = pd.read_csv(f, delimiter=' ').round(5).values.tolist()
        if len(tflist) < 4**12:
            print("Skip %s since it has less rows than 4**12" % pbmname)
            buggedtf += 1
            continue
        for row_key in container:
            seqidx = row_key[3]
            diff = tflist[seqidx][0]
            zscore = tflist[seqidx][1]
            if np.isnan(zscore):
                zscore = 0
            if np.isnan(diff):
                diff = 0
            pval = scipy.stats.norm.sf(abs(zscore))*2
            add = True
            if filteropt == 1:
                # if z-score is chosen then filterval is the maximum of item shown
                if len(container[row_key]) >= filterval:
                    least_idx = min(enumerate(container[row_key]),key=lambda x:abs(x[1][1]))[0]
                    if abs(tflist[seqidx][1]) > abs(container[row_key][least_idx][1]):
                        del container[row_key][least_idx]
                    else:
                        add = False
            # filteropt = 2, if z-score is chosen then filterval is the p-val threshold
            elif pval > filterval:
                    add = False
            # E-score calculation is here
            if add:
                # we don't do escore calculation here
                if spec_ecutoff == -1 or nonspec_ecutoff == -1:
                    container[row_key].append([diff,zscore,pval,"N/A",pbmname])
                else:
                    test_start = timer()
                    # E-score calculation: 0.05 seconds each
                    # For 10k rows, total: 141.34secs, from e-score 128.56331secs
                    # For 50k rows, total: 771.42 secs, from e-score: 752.123secs
                    # another example: 2547.41secs, from e-score: 2523.96897secs
                    eshort_path = "%s/%s_escore.txt" % (config.ESCORE_DIR,pbmname)
                    short2long_map = "%s/index_short_to_long.csv" % (config.ESCORE_DIR)
                    isbound = utils.isbound_escore_18mer(row_key[2], eshort_path, short2long_map, spec_ecutoff, nonspec_ecutoff)
                    container[row_key].append([diff,zscore,pval,isbound,pbmname])
                    test_end = timer()
                    escore_total_time += (test_end-test_start)

        print("Total e-score time %.5f" % escore_total_time)
        ready_count.value += 1
        print("Total running time for {}: {:.2f}secs".format(pbmname,time.time()-start))

    # remove seqidx and 18mer as it is not needed anymore
    newcontainer = {}
    for row_key in container:
        newcontainer[row_key[:-2]] = container[row_key]
    return newcontainer

def format2tbl(tbl,gene_names,filteropt=1):
    '''
    This function saves tbl as csvstring

    Input:
      tbl is a dictionary of (rowidx,seq):[diff,zscore,tfname] or [diff,p-val,escore,tfname]
    '''

    with open(config.PBM_HUGO_MAPPING) as f:
        pbmtohugo = {}
        for line in f:
            linemap = line.strip().split(":")
            pbmtohugo[linemap[0]] = linemap[1].split(",")

    #gapdata = read_gapfile(app.config['GAP_FILE'])

    sorted_key = sorted(tbl.keys())
    datavalues = []
    for row_key in sorted_key:
        if not tbl[row_key]: # probably empty row
            continue
        row = row_key[0]
        seq = row_key[1]
        wild = seq[0:5] + seq[5] + seq[6:11]
        mut = seq[0:5] + seq[11] + seq[6:11]
        sorted_val = sorted(tbl[row_key],reverse=True,key=lambda x:abs(x[1]))
        for row_val in sorted_val: # [diff,zscore,pval,isbound,pbmname]
            rowdict = {'row':row,'wild':wild,'mutant':mut,'diff':row_val[0]}
            pbmname = row_val[4]
            rowdict['z_score'] =  row_val[1]
            rowdict['p_value'] =  row_val[2]
            rowdict['binding_status'] = row_val[3]
            if pbmname  == 'None':
                rowdict['TF_gene'] = ""
                rowdict['pbmname'] = "None"
                #rowdict['gapmodel'] = "None" # vmartin: comment for now
            else:
                rowdict['TF_gene'] = ",".join([gene for gene in pbmtohugo[pbmname] if gene in gene_names])
                rowdict['pbmname'] = pbmname
                #rowdict['gapmodel'] = gapdata[pbmname] # vmartin: comment for now
            datavalues.append(rowdict)

    #colnames = ["row","wild","mutant","diff","z_score","p_value","TF_gene","binding_status","gapmodel","pbmname"]
    colnames = ["row","wild","mutant","diff","z_score","p_value","TF_gene","binding_status","pbmname"]
    return colnames,datavalues

def postprocess(datalist,gene_names,filteropt=1,filterval=1):
    '''
    Aggregate the result from the different processes.
    '''
    maintbl = {}
    for ddict in datalist:
        if not maintbl:
            maintbl = ddict
        else:
            if filteropt == 1: # z-score
                for row_key in ddict:
                    for row_val in ddict[row_key]:
                        least_idx = min(enumerate(maintbl[row_key]),key=lambda x:abs(x[1][1]))[0]
                        # row_val[1] is the t-value
                        if abs(row_val[1]) > abs(maintbl[row_key][least_idx][1]):
                            del maintbl[row_key][least_idx]
                            maintbl[row_key].append(row_val)
            else: # filteropt == 2 -- p-value
                for row_key in ddict:
                    maintbl[row_key].extend(ddict[row_key])
    return format2tbl(maintbl,gene_names,filteropt)

def parse_tfgenes(filepath, prefix = "prediction6mer.", sufix = ".txt"):
    genes = open(filepath).read().splitlines()
    unique_pbms = list({prefix + tf + sufix for gene in genes for tf in config.HUGO_PBM_MAPPING[gene]})
    return {"genes":genes,"pbms":unique_pbms}

def do_prediction(intbl, pbms, gene_names,
                  filteropt="p-value", filterval=0.0001, spec_ecutoff=0.4, nonspec_ecutoff=0.35):
    """
    intbl: preprocessed table
    filteropt: p-value or z-score
    filterval: # TFs for opt z-score and p-val cutoff for p-value
    """

    # intbl: #rowidx,seq,val,diff,t,pbmname,escore_seq
    start_time = time.time()

    # move the comment here for testing
    predfiles = [config.PREDDIR + "/" + pbm for pbm in pbms] # os.listdir(preddir)
    preds = utils.chunkify(predfiles,config.PCOUNT) # chunks the predfiles for each process

    # need to use manager here
    shared_ready_sum = mp.Manager().Value('i', 0)

    pool = mp.Pool(processes=config.PCOUNT)
    if filteropt == "p-value":
        filterval = float(filterval)
    else: #z-score
        filterval = int(filterval)
    async_pools = [pool.apply_async(predict, (preds[i], intbl, shared_ready_sum, filteropt, filterval, spec_ecutoff, nonspec_ecutoff)) for i in range(0,len(preds))]

    total = len(predfiles)
    while not all([p.ready() for p in async_pools]):
        time.sleep(2)

    res = [p.get() for p in async_pools]
    pool.terminate()

    colnames,datavalues = postprocess(res,gene_names,filteropt,filterval)

    return colnames,datavalues

def main():
    """nonspec_bind_cutoff = 0.35
    spec_bind_cutoff = 0.4
    chrver = "hg38"
    filteropt = "p-value"
    filterval = 0.001
    fileinput = "inputs/test.vcf"
    genelist = "inputs/gene_input.txt"
    outpath = "result.csv"""

    parser = argparse.ArgumentParser(description = 'TF Mutation Predictions')
    parser.add_argument('-i', '--inputfile', action="store", dest="inputfile", type=str,
                        help='Input mutation file in .vcf, .tsv, .csv, or .txt format.')
    parser.add_argument('-g', '--genefile', action="store", dest="genefile", type=str,
                        help='A file that contains all TF genes that are desired.')
    parser.add_argument('-t', '--filetype', action="store", dest="filetype", type=str,
                        default="", help='File type can specify: vcf, icgc, customseq, or mut')
    parser.add_argument('-o', '--outpath', action="store", dest="outpath", type=str,
                        default="out.tsv", help='name of the .tsv file that is made')
    parser.add_argument('-f', '--filteropt', action="store", dest="filteropt", type=str,
                        default="p-value", help='p-value or z-score')
    parser.add_argument('-v', '--filterval', action="store", dest="filterval", type=float,
                        default=0.0001, help='# TFs for opt z-score and p-val cutoff for p-value')
    parser.add_argument('-c', '--chrver', action="store", dest="chrver", type=str,
                        default="hg19", help='Chromosome version, can be hg19 or hg38')
    parser.add_argument('-E', '--escorespec', action="store", dest="escorespec", type=float,
                        default=0.4, help='PBM E-score specific cutoff.')
    parser.add_argument('-e', '--escorenonspec', action="store", dest="escorenonspec", type=float,
                        default=0.35, help='PBM E-score non-specific cutoff.')
    args = parser.parse_args()

    #python3 qbic.py -i testing_resources/single_mutations_sample.tsv -g testing_resources/gene_input.txt -c hg19
    # input_mutation_test.vcf
    # -t mut
    #TfX E2F

    if not args.inputfile or not args.genefile:
        raise Exception('-i (--inputfile) and -g (--genefile) are required')

    tbl = inittbl(args.inputfile, args.chrver, filetype = args.filetype)

    input_genes = parse_tfgenes(args.genefile)
    colnames, datavalues = do_prediction(tbl, input_genes["pbms"], input_genes["genes"], args.filteropt, args.filterval,
                                         args.escorespec, args.escorenonspec)

    print("Writing output to %s" % args.outpath)
    pd.DataFrame(datavalues).to_csv(args.outpath, index = False, columns = colnames, sep="\t")

if __name__=="__main__":
    main()
