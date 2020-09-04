import pandas as pd
import time
import os
import numpy as np
import scipy.stats
import argparse
import collections
import sys
import concurrent.futures as cc
import functools as ft
import ctypes # shared arr
from timeit import default_timer as timer

import utils
import config

# comment : df vs dict

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
        df['chromosome'] = df['chromosome'].str.replace("chr", "")

        grouped = df.groupby('chromosome',sort=True)
        dataset = {key:item for key,item in grouped}

        cidx_list = [str(a) for a in list(range(1,23))+['X','Y'] if str(a) in dataset]
        with cc.ThreadPoolExecutor(config.PCOUNT) as executor:
            # vm: first parallel
            result = executor.map(lambda x: chrom_cidx_helper(*x), [(cidx, dataset[cidx], chromosome_version, kmer) for cidx in cidx_list])
        result = [r for group in result for r in group] # flatten list

    result = sorted(result,key=lambda result:result[0])
    # example row in result: [73, 'CCAACCAACCCA', 'ATTCCAACCAACCCCCTA', 5263444, 0, 0, 'None']
    print("Time to preprocess: {:.2f}secs".format(time.time()-start))
    return result

def chrom_cidx_helper(cidx, cidx_dataset, chromosome_version, kmer):
    print("Iterating dataset for chromosome {}...".format(cidx))
    chromosome = utils.get_chrom(config.CHRDIR + "/" + chromosome_version + "/chr." + str(cidx) + '.fa.gz')
    result = []
    for idx,row in cidx_dataset.iterrows():
        pos = row['pos'] - 1
        if row['mutated_from'] != chromosome[pos]:
            error = "For the input mutation %s>%s at position %s in chromosome %s, the mutated_from nucleotide (%s) does not match the nucleotide in the %s reference genome (%s). Please check the input data and verify that the correct version of the reference human genome was used." % (row['mutated_from'], row['mutated_to'], row['pos'], row['chromosome'], row['mutated_from'], chromosome_version, chromosome[pos])
            #raise Exception(error)
            print(error)
        seq = chromosome[pos-kmer+1:pos+kmer] + row['mutated_to'] #-5,+6
        # for escore, just use 8?
        escore_seq = chromosome[pos-9+1:pos+9] + row['mutated_to']
        result.append([idx,seq,escore_seq,utils.seqtoi(seq),0,0,"None"]) #rowidx,seq,escore_seq,val,diff,t,pbmname
    return result

def predict(predlist, dataset, emap,
            filteropt="p-value", filterval=0.001, spec_ecutoff=0.4, nonspec_ecutoff=0.35,
            q=None, num_threads=None):
    """
    for the container list, key is a tuple of: (rowidx,sequence,seqidx)
    and each element in value is a list of: [diff,z-score,pbmname
    If spec_ecutoff and nonspec_ecutoff == -1 then no escore calculation is done

    filteropt = p-value, t-value, none
                if none then just calculate escore

    return:
     filteropt=1: diff,z_score,tfname
     filteropt=2: diff,z_score,p_val,escore,tfname

    update:
     delete ready_count
    """
    #[96, 'TCATGGTGGGTT', GCTTCATGGTGGGTGGAT, 13872815, 0, 0, '-'] -- 37, 'GCCCAGAAAGGA', 9773096
    '''
    EDIT -- this only works for p-values, z-scores need a different memory access pattern
    '''

    # store the results here
    res = []

    # move here so we just use partial function
    pph_partial = ft.partial(pred_helper, **{'dataset':dataset, 'emap':emap, 'filterval':filterval,
                                                'spec_ecutoff':spec_ecutoff, 'nonspec_ecutoff':nonspec_ecutoff})
    # vmartin: confirm how is ready_count used
    if num_threads is  None:
        #iterate for each transcription factor
        for i in range(0,len(predlist)):
            res += [pph_partial(predlist[i])]
            #ready_count.value += 1
    else:
        # concurrent execution for improved I/O,
        # after a partial function is created, mapping is easier
        # not necessarily give benefit but can be
        with cc.ThreadPoolExecutor(max_workers = num_threads) as executor:
            res = executor.map(pph_partial, predlist)
            #ready_count.value += len(predlist)

    # concatenate the individual tf containers
    res = pd.concat(res, axis=0, ignore_index=True)

    # may modify downstream tasks to accept dataframe
    # return res

    # vmartin: are we really using q here?
    if q is None:
        return res
    else:
        q.put(res)
        return None

def pred_helper(pred, dataset, emap, filterval=0.001, filteropt="p-value",spec_ecutoff=0.4, nonspec_ecutoff=0.35):
    '''
    Helper function to shorten predict. Returns a DataFrame of thresholded results
    pred: the 12-mer prediction tables
    '''
    pbmname = '.'.join(pred.split(".")[1:-1])
    print("Processing " + pbmname)
    start = time.time()

    with open(pred, 'r') as f:
        tf_df = pd.read_csv(f, delimiter=' ').round(5).fillna(0)
        if tf_df.shape[0] < 4**12:
            print("Skip %s since it has less rows than 4**12" % pbmname)
            return None

    # create a new local dataframe for computation + filtration
    container = pd.DataFrame(dataset, columns = ['row_key', '12mer', '18mer', 'seqidx', 'diff','z-score', 'pbmname'])

    # extract the z-scores and pvalues
    container['z-score'] = np.array(tf_df.iloc[container['seqidx'], 1]) #.to_numpy() # copy values otherwise pd.Series index issue
    container['p-val'] = scipy.stats.norm.sf(np.abs(container['z-score']))*2

    if filteropt == "p-value":
        # drop the p values above threshold (insignificant)
        container = container[container['p-val'] <= filterval]
    else: # "z-score"
        # if z-score, we take the maximum "filterval"
        # change this to heap
        # nlargest equivalent to below but a bit faster
        container = container.sort_values('z-score',ascending = False).head(filterval)

    # collect diff (done after thresholding)
    container['diff'] = np.array(tf_df.iloc[container['seqidx'], 0]) #.to_numpy() # copy values otherwise pd.Series index issue

    # create a bound column and set the pbmname
    container['binding_status'] = "N/A"
    container['pbmname'] = pbmname

    # resort columns
    container = container[['row_key', '12mer', '18mer', 'seqidx', 'diff','z-score', 'p-val', 'binding_status', 'pbmname']]

    if spec_ecutoff != -1 and nonspec_ecutoff != -1:
        eshort_path = "%s/%s_escore.txt" % (config.ESCORE_DIR,pbmname)
        eshort = pd.read_csv(eshort_path, header=None, index_col=None, dtype=np.float32)[0] # pd.Series -- remove from utils
        elong = np.array(eshort.iloc[emap]) #.to_numpy() # Moved out of isbound_escore_18mer()
        container['binding_status'] = container['18mer'].apply(lambda x: utils.isbound_escore_18mer(x, elong, spec_ecutoff, nonspec_ecutoff))

    container.drop(columns=['seqidx', '18mer'], inplace=True)

    print("Total running time for {}: {:.2f}secs".format(pbmname,time.time()-start))

    return container

def postprocess(datalist,predfiles,gene_names,filteropt=1,filterval=1):
    '''
    Aggregate the result from the different processes.

    TODO -- z-score
    '''
    datalist = pd.concat(datalist, ignore_index=True, axis=0)

    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)

    if filteropt == "p-value":
        datalist.sort_values(by = ['row_key', '12mer', 'p-val'], ascending=True, inplace=True)
    else: # if z-score then we need to reduce here
        datalist = datalist.iloc[(-datalist['z-score'].abs()).argsort()].groupby('row_key').head(filterval)
        datalist.sort_values(by = ['row_key', '12mer'], ascending=True, inplace=True)

    # read in pbm to hugo map and split
    pbmtohugo = pd.read_csv(config.PBM_HUGO_MAPPING, sep=':', index_col=0, header=None)[1].str.split(',')

    # reconstruct wild type and mutant strings
    datalist['wild'] = datalist['12mer'].str.slice(stop=11)
    datalist['mutant'] = datalist['12mer'].str.slice(stop=5) + datalist['12mer'].str.get(11)+ datalist['12mer'].str.slice(start=6, stop=11)

    gene_names_set = set(gene_names)
    predfiles = ['.'.join(p.split(".")[1:-1]) for p in predfiles]
    tf_gene_dict = {p: ",".join([gene for gene in pbmtohugo[p] if gene in gene_names_set]) for p in predfiles}
    datalist['TF_gene'] = datalist['pbmname'].apply(lambda x: tf_gene_dict.get(x, x))

    # reindex and rename the columns
    datalist = datalist[['row_key', 'wild', 'mutant', 'diff', 'z-score', 'p-val', 'TF_gene', 'binding_status', 'pbmname']]
    datalist.columns = ["row","wild","mutant","diff","z_score","p-value","TF_gene","binding_status","pbmname"]

    return datalist

def parse_tfgenes(filepath, prefix = "prediction6mer.", sufix = ".txt"):
    genes = open(filepath).read().splitlines()
    unique_pbms = list({prefix + tf + sufix for gene in genes for tf in config.HUGO_PBM_MAPPING[gene]})
    return {"genes":genes,"pbms":unique_pbms}

def do_prediction(intbl, pbms, gene_names,
                  filteropt="p-value", filterval=0.0001, spec_ecutoff=0.4, nonspec_ecutoff=0.35, num_threads=None):
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

    if filteropt == "p-value":
        filterval = float(filterval)
    else: #z-score
        filterval = int(filterval)

    # collect the short2long_map -- shared, so only one i/o
    emap = pd.read_csv("%s/index_short_to_long.csv" % (config.ESCORE_DIR), header=0, index_col=0, sep=',', dtype='Int32') # pd.DataFrame
    emap = np.array(emap[emap.columns[0]]) - 1 #emap[emap.columns[0]].to_numpy() - 1


    # --- PARALLEL PART ---
    # need to use manager here
    #shared_ready_sum = mp.Manager().Value('i', 0)
    # prepare all parameters but predlist
    predict_partial = ft.partial(predict, **{'dataset':intbl, 'emap':emap,
            'filteropt':filteropt, 'filterval':filterval, 'spec_ecutoff':spec_ecutoff, 'nonspec_ecutoff':nonspec_ecutoff, 'num_threads':num_threads})

    with cc.ProcessPoolExecutor(config.PCOUNT) as executor:
        res = executor.map(predict_partial, preds)

    return postprocess(res,predfiles,gene_names,filteropt,filterval)

def main():
    """nonspec_bind_cutoff = 0.35
    spec_bind_cutoff = 0.4
    chrver = "hg38"
    filteropt = "p-value"
    filterval = 0.001
    fileinput = "inputs/test.vcf"
    genelist = "inputs/gene_input.txt"
    outpath = "result.csv"""

    # cProfile capitulation
    # https://stackoverflow.com/questions/53890693/cprofile-causes-pickling-error-when-running-multiprocessing-python-code
    import cProfile
    if sys.modules['__main__'].__file__ == cProfile.__file__:
        print("here cprofile")
        import qbic # re-import main (does *not* use cache or execute as __main__)
        globals().update(vars(qbic))  # Replaces current contents with newly imported stuff
        sys.modules['__main__'] = qbic  # Ensures pickle lookups on __main__ find matching version


    parser = argparse.ArgumentParser(description = 'TF Mutation Predictions')
    parser.add_argument('-i', '--inputfile', action="store", dest="inputfile", type=str,
                        help='Input mutation file in .vcf, .tsv, .csv, or .txt format.')
    parser.add_argument('-g', '--genesfile', action="store", dest="genesfile", type=str,
                        help='A file that contains all TF genes that are desired.')
    parser.add_argument('-t', '--filetype', action="store", dest="filetype", type=str,
                        default="", help='File type can specify: vcf, icgc, customseq, or mut')
    parser.add_argument('-o', '--outpath', action="store", dest="outpath", type=str,
                        default="out.tsv", help='name of the .tsv file that is made')
    parser.add_argument('-f', '--filteropt', action="store", dest="filteropt", type=str,
                        default="p-value", choices=['p-value', 'z-score'], help='p-value or z-score')
    parser.add_argument('-v', '--filterval', action="store", dest="filterval", type=float,
                        default=0.0001, help='# TFs for opt z-score and p-val cutoff for p-value')
    parser.add_argument('-c', '--chrver', action="store", dest="chrver", type=str,
                        default="hg19", help='Chromosome version, can be hg19 or hg38')
    parser.add_argument('-E', '--escorespec', action="store", dest="escorespec", type=float,
                        default=0.4, help='PBM E-score specific cutoff.')
    parser.add_argument('-e', '--escorenonspec', action="store", dest="escorenonspec", type=float,
                        default=0.35, help='PBM E-score non-specific cutoff.')
    parser.add_argument('-n', '--numthreads', action="store", dest="numthreads", type=int,
                        default=None, help='Number of concurrent file I/O threads to use (per core)')
    args = parser.parse_args()

    #python3 qbic.py -i testing_resources/input_mutation_test.vcf -g testing_resources/gene_input.txt -c hg19

    # vm: added checking
    if not args.inputfile or not args.genesfile:
        raise Exception('-i (--inputfile) and -g (--genefile) are required')
    if args.filterval == "p-value" and (float(args.filterval) < 0 or float(args.filterval) > 1):
        raise Exception('p-value should be between 0 and 1')
    if args.filterval == "z-score" and int(args.filterval) < 1:
        raise Exception('please specify integer for filterval')

    tbl = inittbl(args.inputfile, args.chrver, filetype = args.filetype)
    input_genes = parse_tfgenes(args.genesfile)
    res = do_prediction(tbl, input_genes["pbms"], input_genes["genes"], args.filteropt, args.filterval,
                                         args.escorespec, args.escorenonspec, args.numthreads)
    print("Writing output to %s" % args.outpath)
    res.to_csv(args.outpath, index = False, sep="\t")

if __name__=="__main__":
    main()
