
from celery import Celery,chain
from app import app,celery,mongodb,redisdb
import pymongo

import redisearch
import time,sys
import pandas as pd
import billiard as mp #multiprocessing substitute to enable daemon
import scipy.stats
import numpy as np
import functools as ft

import os

sys.path.insert(0, 'app')
import controller.utils as utils

from timeit import default_timer as timer

import concurrent.futures as cc

# input: table
@celery.task(bind=True)
def inittbl(self,filename,cpath):
    self.update_state(state='PROGRESS',
                  meta={'current': 0, 'total': 1, 'status': 'Preprocessing input...'})
    kmer = 6
    start = time.time()
    file_extension = os.path.splitext(filename)[1]

    result = []
    error = ""

    # TODO: if fast enough, we can also put error checking in here
    if file_extension == ".txt":
        with open(filename) as f:
            idx = 0
            for line in f:
                if "\t" in line:
                    line = line.split("\t")
                else:
                    line = line.split()
                idx += 1
                # line[1] is the base mid nucleotide mutated to
                escore_seq = line[0] + line[1]
                mid_seq = escore_seq[len(escore_seq)//2-6:len(escore_seq)//2+5] + line[1] # the 12mer seq
                result.append([idx,mid_seq,escore_seq,utils.seqtoi(mid_seq),0,0,"None"])
    else:
        if file_extension == ".vcf":
            df = pd.read_csv(filename,sep="\t",header=None).drop(2,1)
            df = df.rename(columns={0:"chromosome",1:"pos",3:"mutated_from",4:"mutated_to"})
            df['chromosome'] = df['chromosome'].map(lambda x:x.replace("chr",""))
        else:
            if file_extension == ".tsv":
                separator = "\t"
            else: # must be csv since we checked it, TODO: can also return error here
                separator = ","
            df = pd.read_csv(filename, sep=separator)
            # if icgc then only take a subset of the columns
            if set(['chromosome','chromosome_start','mutation_type','mutated_from_allele','mutated_to_allele']).issubset(df.columns):
                df = df[['chromosome','chromosome_start','mutation_type','mutated_from_allele','mutated_to_allele']]
                df = df[df['mutation_type'].apply(lambda x: "single base substitution" == x)].drop('mutation_type',1).drop_duplicates() # only take single base mutation
                df = df.rename(columns={"chromosome_start":"pos","mutated_from_allele":"mutated_from","mutated_to_allele":"mutated_to"})
            else: # ['chromosome', 'chromosome_pos', 'mutated_from', 'mutated_to']
                df = df.rename(columns={"chromosome_pos":"pos","mutated_from_allele":"mutated_from","mutated_to_allele":"mutated_to"})
        grouped = df.groupby('chromosome',sort=True)
        dataset = {str(key):item for key,item in grouped}

        for cidx in [str(a) for a in range(1,23)] + ['X','Y']:
            self.update_state(state='PROGRESS',
                      meta={'current': 0, 'total': 1, 'status': 'Preprocessing input for chromosome {}...'.format(cidx)})
            if cidx not in dataset:
                continue
            print("Iterating dataset for chromosome {}...".format(cidx))
            chromosome = utils.get_chrom(cpath + "/chr." + str(cidx) + '.fa.gz')
            for idx,row in dataset[cidx].iterrows():
                pos = row['pos'] - 1
                if row['mutated_from'] != chromosome[pos]:
                    cver = cpath.split("/")[-1]
                    error = "For the input mutation %s>%s at position %s in chromosome %s, the mutated_from nucleotide (%s) does not match the nucleotide in the %s reference genome (%s). Please check the input data and verify that the correct version of the reference human genome was selected in the Data Submission Form." % (row['mutated_from'], row['mutated_to'], row['pos'], row['chromosome'], row['mutated_from'], cver, chromosome[pos])
                    #error = "Found mismatch in the mutation: chromosome %s pos %s mutated_from: %s; but expected: %s. Input mutation coordinate is probably incorrect or different genome version is probably used.\n" % (row['chromosome'],row['pos'],row['mutated_from'],chromosome[pos])
                    break
                seq = chromosome[pos-kmer+1:pos+kmer] + row['mutated_to'] #-5,+6
                # for escore, just use 8?
                escore_seq = chromosome[pos-9+1:pos+9] + row['mutated_to']
                result.append([idx,seq,escore_seq,utils.seqtoi(seq),0,0,"None"]) #rowidx,seq,escore_seq,val,diff,t,pbmname
            if error:
                break
        # cidx_list = [str(a) for a in range(1,23)] + ['X','Y']
        # with cc.ThreadPoolExecutor(app.config["PCOUNT"]) as executor:
        #     # vm: first parallel
        #     result = executor.map(lambda x: chrom_cidx_helper(*x), [(cidx, dataset[cidx], cpath, kmer) for cidx in cidx_list if cidx in dataset])
        # result = sorted(result,key=lambda result:result[0])

    # finish parsing the file, delete it
    if filename.startswith(app.config['UPLOAD_FOLDER']):
        utils.delete_file(filename)

    if error:
        return error
    else:
        result = sorted(result,key=lambda result:result[0])
        # example row in result: [73, 'CCAACCAACCCA', 'ATTCCAACCAACCCCCTA', 5263444, 0, 0, 'None']
        print("Time to preprocess: {:.2f}secs".format(time.time()-start))
        return result

# still doesn't work, TODO
@celery.task(bind=True)
def chrom_cidx_helper(self, cidx, cidx_dataset, chromosome_path, kmer):
    self.update_state(state='PROGRESS',
              meta={'current': 0, 'total': 1, 'status': 'Preprocessing input for chromosome {}...'.format(int(cidx))})
    print("Iterating dataset for chromosome {}...".format(cidx))
    chromosome = utils.get_chrom(chromosome_path + "/chr." + str(cidx) + '.fa.gz')
    result = []
    error = ""
    for idx,row in cidx_dataset.iterrows():
        pos = row['pos'] - 1
        if row['mutated_from'] != chromosome[pos]:
            cver = cpath.split("/")[-1]
            error = "For the input mutation %s>%s at position %s in chromosome %s, the mutated_from nucleotide (%s) does not match the nucleotide in the %s reference genome (%s). Please check the input data and verify that the correct version of the reference human genome was used." % (row['mutated_from'], row['mutated_to'], row['pos'], row['chromosome'], row['mutated_from'], cver, chromosome[pos])
        seq = chromosome[pos-kmer+1:pos+kmer] + row['mutated_to'] #-5,+6
        # for escore, just use 8?
        escore_seq = chromosome[pos-9+1:pos+9] + row['mutated_to']
        result.append([idx,seq,escore_seq,utils.seqtoi(seq),0,0,"None"]) #rowidx,seq,escore_seq,val,diff,t,pbmname
    if error:
        return error
    else:
        return result

#==================================== Prediction Part ====================================

def predict(predlist, dataset, ready_count, emap,
            filteropt="p-value", filterval=0.001, spec_ecutoff=0.4, nonspec_ecutoff=0.35):
    """
    for the container list, key is a tuple of: (rowidx,sequence,seqidx)
    and each element in value is a list of: [diff,z-score,pbmname]
    If spec_ecutoff and nonspec_ecutoff == -1 then no escore calculation is done

    filteropt = p-value, t-value, none
                if none then just calculate escore

    return:
     filteropt=1: diff,z_score,tfname
     filteropt=2: diff,z_score,p_val,escore,tfname
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
    #iterate for each transcription factor
    for i in range(0,len(predlist)):
        res += [pph_partial(predlist[i])]
        ready_count.value += 1

    # concatenate the individual tf containers
    res = pd.concat(res, axis=0, ignore_index=True)

    # may modify downstream tasks to accept dataframe
    # return res
    return res

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
        eshort_path = "%s/%s_escore.txt" % (app.config["ESCORE_DIR"],pbmname)
        eshort = pd.read_csv(eshort_path, header=None, index_col=None, dtype=np.float32)[0] # pd.Series -- remove from utils
        elong = np.array(eshort.iloc[emap]) #.to_numpy() # Moved out of isbound_escore_18mer()
        container['binding_status'] = container['18mer'].apply(lambda x: utils.isbound_escore_18mer(x, elong, spec_ecutoff, nonspec_ecutoff))

    container.drop(columns=['seqidx', '18mer'], inplace=True)

    print("Total running time for {}: {:.2f}secs".format(pbmname,time.time()-start))

    return container

def read_gapfile(gapfile):
    df = pd.read_csv(gapfile)
    return dict(zip(df.upbm_filenames, df.gapmodel))
#
# def format2tbl(tbl,gene_names,filteropt=1):
#     '''
#     This function saves tbl as csvstring
#
#     Input:
#       tbl is a dictionary of (rowidx,seq):[diff,zscore,tfname] or [diff,p-val,escore,tfname]
#     '''
#
#     with open(app.config['PBM_HUGO_MAPPING']) as f:
#         pbmtohugo = {}
#         for line in f:
#             linemap = line.strip().split(":")
#             pbmtohugo[linemap[0]] = linemap[1].split(",")
#
#     #gapdata = read_gapfile(app.config['GAP_FILE'])
#
#     sorted_key = sorted(tbl.keys())
#     datavalues = []
#     for row_key in sorted_key:
#         if not tbl[row_key]: # probably empty row
#             continue
#         row = row_key[0]
#         seq = row_key[1]
#         wild = seq[0:5] + seq[5] + seq[6:11]
#         mut = seq[0:5] + seq[11] + seq[6:11]
#         sorted_val = sorted(tbl[row_key],reverse=True,key=lambda x:abs(x[1]))
#         for row_val in sorted_val: # [diff,zscore,pval,isbound,pbmname]
#             rowdict = {'row':row,'wild':wild,'mutant':mut,'diff':row_val[0]}
#             pbmname = row_val[4]
#             rowdict['z_score'] =  row_val[1]
#             rowdict['p_value'] =  row_val[2]
#             rowdict['binding_status'] = row_val[3]
#             if pbmname  == 'None':
#                 rowdict['TF_gene'] = ""
#                 rowdict['pbmname'] = "None"
#                 #rowdict['gapmodel'] = "None" # vmartin: comment for now
#             else:
#                 rowdict['TF_gene'] = ",".join([gene for gene in pbmtohugo[pbmname] if gene in gene_names])
#                 rowdict['pbmname'] = pbmname
#                 #rowdict['gapmodel'] = gapdata[pbmname] # vmartin: comment for now
#             datavalues.append(rowdict)
#
#     #colnames = ["row","wild","mutant","diff","z_score","p_value","TF_gene","binding_status","gapmodel","pbmname"]
#     colnames = ["row","wild","mutant","diff","z_score","p_value","TF_gene","binding_status","pbmname"]
#     return colnames,datavalues
#
# def postprocess(datalist,gene_names,filteropt=1,filterval=1):
#     '''
#     Aggregate the result from the different processes.
#     '''
#     maintbl = {}
#     for ddict in datalist:
#         if not maintbl:
#             maintbl = ddict
#         else:
#             if filteropt == 1: # z-score
#                 for row_key in ddict:
#                     for row_val in ddict[row_key]:
#                         least_idx = min(enumerate(maintbl[row_key]),key=lambda x:abs(x[1][1]))[0]
#                         # row_val[1] is the t-value
#                         if abs(row_val[1]) > abs(maintbl[row_key][least_idx][1]):
#                             del maintbl[row_key][least_idx]
#                             maintbl[row_key].append(row_val)
#             else: # filteropt == 2 -- p-value
#                 for row_key in ddict:
#                     maintbl[row_key].extend(ddict[row_key])
#     return format2tbl(maintbl,gene_names,filteropt)

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
    pbmtohugo = pd.read_csv(app.config["PBM_HUGO_MAPPING"], sep=':', index_col=0, header=None)[1].str.split(',')

    # reconstruct wild type and mutant strings
    datalist['wild'] = datalist['12mer'].str.slice(stop=11)
    datalist['mutant'] = datalist['12mer'].str.slice(stop=5) + datalist['12mer'].str.get(11)+ datalist['12mer'].str.slice(start=6, stop=11)

    gene_names_set = set(gene_names)
    predfiles = ['.'.join(p.split(".")[1:-1]) for p in predfiles]
    tf_gene_dict = {p: ",".join([gene for gene in pbmtohugo[p] if gene in gene_names_set]) for p in predfiles}
    datalist['TF_gene'] = datalist['pbmname'].apply(lambda x: tf_gene_dict.get(x, x))

    # reindex and rename the columns
    datalist = datalist[['row_key', 'wild', 'mutant', 'diff', 'z-score', 'p-val', 'TF_gene', 'binding_status', 'pbmname']]
    datalist.columns = ["row","wild","mutant","diff","z_score","p_value","TF_gene","binding_status","pbmname"]

    return datalist

#==========================================================
@celery.task()
def drop_collection(task_id):
    '''
    Make this a celery task so we can schedule it -- done?
    '''
    print("Remove key/index for %s from db"%task_id)
    """
    client = redisearch.Client(task_id)
    client.drop_index()
    db.delete(task_id)
    db.delete("%s:cols"%task_id)
    """
    collection = mongodb[task_id]
    collection.drop()

def savetomongo(req_id,datavalues,expired_time):
    """
    datavalues: list of dictionary of values to be put in the database
    """
    collection = mongodb[req_id]
    #collection.insert_many(df.to_dict("records"))
    indexes = []
    collection.insert_many(datavalues)
    collection.create_index([
                    ("row", pymongo.ASCENDING), ("row", pymongo.DESCENDING),
                    ("diff", pymongo.ASCENDING), ("diff", pymongo.DESCENDING),
                    ("z_score", pymongo.ASCENDING), ("z_score", pymongo.DESCENDING),
                    ("p_value", pymongo.ASCENDING), ("p_value", pymongo.DESCENDING),
                    ("wild", pymongo.TEXT), ("mutant", pymongo.TEXT)
                ])
    print("Save collection %s to database" % req_id)
    # ---- set expiry for columns and documents ----
    drop_collection.apply_async((req_id,), countdown=expired_time)

#https://github.com/MehmetKaplan/Redis_Table
@celery.task(bind=True)
def do_prediction(self, intbl, selections, gene_names,
                  filteropt=1, filterval=1, spec_ecutoff=0.4, nonspec_ecutoff=0.35):
    '''
    intbl: preprocessed table
    filteropt: 1 for highest t-val, 2 for p-val cutoff
    filterval: # TFs for opt 1 and p-val cutoff for opt 2
    '''

    if type(intbl) is str: # got an error in the pipeline from inittbl
        return {'current': 1, 'total': 1, 'error': intbl}

    # intbl: #rowidx,seq,val,diff,t,pbmname,escore_seq
    start_time = time.time()

    #while not inittask.ready():
    #    time.sleep(1)
    #intbl = inittask.get()
    predfiles = [app.config['PREDDIR'] + "/" + s for s in selections] # os.listdir(preddir)
    preds = [l for l in utils.chunkify(predfiles,app.config['PCOUNT']) if len(l) != 0] # chunks the predfiles for each process

    # collect the short2long_map -- shared, so only one i/o
    emap = pd.read_csv("%s/index_short_to_long.csv" % (app.config["ESCORE_DIR"]), header=0, index_col=0, sep=',', dtype='Int32') # pd.DataFrame
    emap = np.array(emap[emap.columns[0]]) - 1 #emap[emap.columns[0]].to_numpy() - 1

    # ---- MULTIPROCESSING PART ----
    pool = mp.Pool(processes=app.config['PCOUNT'])
    # need to use manager here
    shared_ready_sum = mp.Manager().Value('i', 0)

    predict_partial = ft.partial(predict, **{'dataset':intbl, 'ready_count':shared_ready_sum, 'emap':emap,
            'filteropt':filteropt, 'filterval':filterval, 'spec_ecutoff':spec_ecutoff, 'nonspec_ecutoff':nonspec_ecutoff})
    async_pools = [pool.apply_async(predict_partial, (preds[i], )) for i in range(0,len(preds))]

    # run the job, update progress bar
    total = len(predfiles)
    while not all([p.ready() for p in async_pools]):
        time.sleep(2) # super important to avoid checking every loop
        self.update_state(state='PROGRESS',
                          meta={'current': shared_ready_sum.value, 'total': total, 'status': 'Processing input data...'})

    res = [p.get() for p in async_pools]
    self.update_state(state='PROGRESS',
                          meta={'current': shared_ready_sum.value, 'total': total, 'status': 'post-processing'})
    print("Terminate all children process..")
    pool.terminate() # terminate to kill all child processes !!! Like.. super important,
                     # to avoid memory leak, seriously...
    datavalues = postprocess(res,predfiles,gene_names,filteropt,filterval)

    ''' SET the values in redis '''
    #print("marktesting",colnames,datavalues)
    savetomongo(self.request.id, datavalues.to_dict('records') ,app.config['USER_DATA_EXPIRY'])
    # significance_score can be z-score or p-value depending on the out_type

    #db.expire("%s:vals:*" % self.request.id, app.config['USER_DATA_EXPIRY'])

    return {'current': shared_ready_sum.value, 'total': len(predfiles), 'status': 'Task completed!',
            'result': 'done', 'taskid': self.request.id,
            'time':(time.time()-start_time)} # -- somehow cannot do jsonify(postproc)
