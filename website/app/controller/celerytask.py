
from celery import Celery,chain
from app import app,celery,db

import redisearch
import time,sys
import pandas as pd
import billiard as mp #multiprocessing substitute to enable daemon
import scipy.stats
import numpy as np

import os

sys.path.insert(0, 'app')
import controller.utils as utils

from timeit import default_timer as timer

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
            df = pd.read_csv(filename,
                        sep=separator,
                        usecols=['chromosome','chromosome_start','mutation_type','mutated_from_allele','mutated_to_allele'])
            df = df[df['mutation_type'].apply(lambda x: "single base substitution" == x)].drop('mutation_type',1).drop_duplicates() # only take single base mutation
            df = df.rename(columns={"chromosome_start":"pos","mutated_from_allele":"mutated_from","mutated_to_allele":"mutated_to"})

        grouped = df.groupby('chromosome',sort=True)
        dataset = {key:item for key,item in grouped}

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
                    error = "Found mismatch in the mutation: chromosome %s pos %s mutated_from: %s; but expected: %s. Input mutation coordinate is probably incorrect or different genome version is probably used.\n" % (row['chromosome'],row['pos'],row['mutated_from'],chromosome[pos])
                    break
                seq = chromosome[pos-kmer+1:pos+kmer] + row['mutated_to'] #-5,+6
                # for escore, just use 8?
                escore_seq = chromosome[pos-9+1:pos+9] + row['mutated_to']
                result.append([idx,seq,escore_seq,utils.seqtoi(seq),0,0,"None"]) #rowidx,seq,escore_seq,val,diff,t,pbmname
            if error:
                break

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

#==================================== Prediction Part ====================================

def predict(predlist, dataset, ready_count,
            filteropt=1, filterval=1, spec_ecutoff=0.4, nonspec_ecutoff=0.35):
    '''
    for the container list, key is a tuple of: (rowidx,sequence,seqidx)
    and each element in value is a list if: [diff,z-score,pbmname]

    return:
     filteropt=1: diff,z_score,tfname
     filteropt=2: diff,p_val,escore,tfname
    '''
    buggedtf = 0
    #[96, 'TCATGGTGGGTT', GCTTCATGGTGGGTGGAT, 13872815, 0, 0, '-'] -- 37, 'GCCCAGAAAGGA', 9773096
    if filteropt == 1: #t-value
        container = {tuple(row[:4]):[[0,0,1,"None","None"]] for row in dataset} # rowidx,12mer,18mer,seqidx : [diff,z,p,bind,pbmname]
    else: #p-value
        # leave this empty as for p-value, we don't have to compare and the size is dynamic
        container = {tuple(row[:4]):[] for row in dataset}

    test_total_time = 0
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
                if spec_ecutoff == -1 or nonspec_ecutoff == -1:
                    container[row_key].append([diff,zscore,pval,"N/A",pbmname])
                else:
                    test_start = timer()
                    # E-score calculation: 0.05 seconds each
                    # For 10k rows, total: 141.34secs, from e-score 128.56331secs
                    # For 50k rows, total: 771.42 secs, from e-score: 752.123secs
                    # another example: 2547.41secs, from e-score: 2523.96897secs
                    isbound = utils.isbound_escore_18mer(row_key[2],pbmname,app.config['ESCORE_DIR'],spec_ecutoff,nonspec_ecutoff)
                    container[row_key].append([diff,zscore,pval,isbound,pbmname])
                    test_end = timer()
                    test_total_time += (test_end-test_start)

        print("Total e-score time %.5f" % test_total_time)
        ready_count.value += 1
        print("Total running time for {}: {:.2f}secs".format(pbmname,time.time()-start))

    # remove seqidx and 18mer as it is not needed anymore
    newcontainer = {}
    for row_key in container:
        newcontainer[row_key[:-2]] = container[row_key]
    return newcontainer

def read_gapfile(gapfile):
    df = pd.read_csv(gapfile)
    return dict(zip(df.upbm_filenames, df.gapmodel))

def format2tbl(tbl,gene_names,filteropt=1):
    '''
    This function saves tbl as csvstring

    Input:
      tbl is a dictionary of (rowidx,seq):[diff,zscore,tfname] or [diff,p-val,escore,tfname]
    '''

    with open(app.config['PBM_HUGO_MAPPING']) as f:
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

#==========================================================
@celery.task()
def drop_index(task_id):
    '''
    Make this a celery task so we can schedule it
    '''
    print("Remove key/index for %s from redis"%task_id)
    client = redisearch.Client(task_id)
    client.drop_index()
    db.delete(task_id)
    db.delete("%s:cols"%task_id)

def savetoredis(req_id,colnames,datavalues,expired_time):
    db.hmset("%s:cols"%req_id,{'cols':colnames})
    client = redisearch.Client(req_id)
    indexes = []
    for col in colnames:
        if "score" in col or "diff" in col or "row" in col or "z_score" in col or "p_value" in col:
            indexes.append(redisearch.NumericField(col,sortable=True))
        else:
            indexes.append(redisearch.TextField(col,sortable=True))
    client.create_index(indexes)
    for i in range(0,len(datavalues)):
        fields = {colnames[j]:datavalues[i][colnames[j]] for j in range(0,len(colnames))}
        client.add_document("%s_%d"%(req_id,i), **fields)
    # ---- set expiry for columns and documents ----
    #db.expire("%s:cols"%req_id,expired_time) let's comment for now and see how it goes
    drop_index.apply_async((req_id,), countdown=expired_time)

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

    # move the comment here for testing
    pool = mp.Pool(processes=app.config['PCOUNT'])
    predfiles = [app.config['PREDDIR'] + "/" + s for s in selections] # os.listdir(preddir)
    preds = utils.chunkify(predfiles,app.config['PCOUNT']) # chunks the predfiles for each process

    # need to use manager here
    shared_ready_sum = mp.Manager().Value('i', 0)

    async_pools = [pool.apply_async(predict, (preds[i], intbl, shared_ready_sum, filteropt, filterval, spec_ecutoff, nonspec_ecutoff)) for i in range(0,len(preds))]

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
    colnames,datavalues = postprocess(res,gene_names,filteropt,filterval)

    ''' SET the values in redis '''
    savetoredis(self.request.id,colnames,datavalues,app.config['USER_DATA_EXPIRY'])
    # significance_score can be z-score or p-value depending on the out_type

    #db.expire("%s:vals:*" % self.request.id, app.config['USER_DATA_EXPIRY'])

    return {'current': shared_ready_sum.value, 'total': len(predfiles), 'status': 'Task completed!',
            'result': 'done', 'taskid': self.request.id,
            'time':(time.time()-start_time)} # -- somehow cannot do jsonify(postproc)
