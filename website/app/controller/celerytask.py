
from celery import Celery,chain
from app import app,celery,db

from decimal import Decimal

import time,sys
import pandas as pd
import billiard as mp #multiprocessing substitute to enable daemon
import scipy.stats

import os

sys.path.insert(0, 'app')
import controller.utils as utils

# input: table
@celery.task(bind=True)
def inittbl(self,filename,cpath):
    self.update_state(state='PROGRESS',
                  meta={'current': 0, 'total': 1, 'status': 'Preprocessing input...'})
    kmer = 6
    start = time.time()
    file_extension = os.path.splitext(filename)[1]
    if file_extension == ".tsv":
        separator = "\t"
    else: # must be csv since we checked it
        separator = ","
    tsv = pd.read_csv(filename,
                sep=separator,
                usecols=['chromosome','chromosome_start','mutation_type','mutated_from_allele','mutated_to_allele'])
    tsv = tsv[tsv['mutation_type'].apply(lambda x: "single base substitution" == x)].drop('mutation_type',1).drop_duplicates() # only take single base mutation
    grouped = tsv.groupby('chromosome',sort=True)
    dataset = {key:item for key,item in grouped}

    result = list()
    error = ""
    for cidx in [str(a) for a in range(1,23)] + ['X','Y']:
        self.update_state(state='PROGRESS',
                  meta={'current': 0, 'total': 1, 'status': 'Preprocessing input for chromosome {}...'.format(cidx)})
        if cidx not in dataset:
            continue
        print("Iterating dataset for chromosome {}...".format(cidx))
        chromosome = utils.get_chrom(cpath + "/chr." + str(cidx) + '.fa.gz')
        for idx,row in dataset[cidx].iterrows():
            pos = row['chromosome_start'] - 1
            if row['mutated_from_allele'] != chromosome[pos]:
                error = "Found mismatch in the mutation: \n{}".format(row)
                break
            seq = chromosome[pos-kmer+1:pos+kmer] + row['mutated_to_allele'] #-5,+6
            # for escore, just use 8?
            esccore_seq = chromosome[pos-9+1:pos+9] + row['mutated_to_allele']
            result.append([idx,seq,esccore_seq,utils.seqtoi(seq),0,0,"-"]) #rowidx,seq,escore_seq,val,diff,t,pbmname
        if error:
            break

    # we finish parsing the file, delete it
    delete_file(filename)
    if error:
        return error
    else:
        result = sorted(result,key=lambda result:result[0])
        print("Time to preprocess: {:.2f}secs".format(time.time()-start))
        return result

#==================================== Prediction Part ====================================

# TODO: remove sharedlist if not needed anymore
def predict(predlist,dataset,sharedlist,filteropt=1,filterval=1):
    '''
    for the container list, key is a tuple of: (rowidx,sequence,seqidx)
    and each element in value is a list if: [diff,z-score,pbmname]

    return:
     filteropt=1: diff,zscore,tfname
     filteropt=2: diff,p-val,escore,tfname
    '''
    buggedtf = 0
    #[96, 'TCATGGTGGGTT', GCTTCATGGTGGGTGGAT, 13872815, 0, 0, '-'] -- 37, 'GCCCAGAAAGGA', 9773096
    if filteropt == 1: #t-value
        container = {tuple(row[:4]):[[0,0,"None"]] for row in dataset} # rowidx,12mer,18mer,seqidx : 0,0,"-"
    else: #p-value
        # leave this empty as for p-value, we don't have to compare and the size is dynamic
        container = {tuple(row[:4]):[] for row in dataset}

    # iterate for each transcription factor
    for i in range(0,len(predlist)):
        start = time.time()
        pbmname = '.'.join(map(str,predlist[i].split(".")[1:-1]))
        print("Processing ",pbmname)
        with open(predlist[i], 'r') as f:
            tflist = pd.read_csv(f, delimiter=' ').round(5).values.tolist()
        if len(tflist) < 4**12:
            print("Skip %s since it has less rows than 4**12" % pbmname)
            buggedtf += 1
            continue
        for row_key in container:
            seqidx = row_key[3]
            if filteropt == 1:
                # if z-score is chosen then filterval is the maximum of item shown
                if len(container[row_key]) >= filterval:
                    least_idx = min(enumerate(container[row_key]),key=lambda x:abs(x[1][1]))[0]
                    if abs(tflist[seqidx][1]) > abs(container[row_key][least_idx][1]):
                        del container[row_key][least_idx]
                        container[row_key].append(tflist[seqidx] + [pbmname])
                else:
                    container[row_key].append(tflist[seqidx] + [pbmname])
            else: # filteropt == 2
                zscore = tflist[seqidx][1]
                pval = scipy.stats.norm.sf(abs(zscore))*2
                # if z-score is chosen then filterval is the p-val threshold
                if pval <= filterval:
                    isbound = utils.isbound_escore_18mer(row_key[2],pbmname,app.config['ESCORE_DIR'])
                    # tflist[seqidx][:-1] -> just diff
                    container[row_key].append(tflist[seqidx][:-1] + [pval,isbound,pbmname])
        sharedlist.append(pbmname) # TODO: delete this
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

    if filteropt == 1:
        colnames = ["rowidx","wild,mutant","diff","z-score","pbmname","TF_gene"]
        metrics = 'z-score'
    else: #filteropt == 1:
        colnames = ["rowidx","wild,mutant","diff","p-value","binding_status","pbmname","TF_gene"]
        metrics = 'p-value'

    with open(app.config['PBM_HUGO_MAPPING']) as f:
        pbmtohugo = {}
        for line in f:
            linemap = line.strip().split(":")
            pbmtohugo[linemap[0]] = linemap[1].split(",")

    sorted_key = sorted(tbl.keys())
    datavalues = []
    for row_key in sorted_key:
        if not tbl[row_key]: # probably empty row
            continue
        seq = row_key[1]
        wild = seq[0:5] + '<span class="bolded-red">' + seq[5] + '</span>' + seq[6:11]
        mut = seq[0:5] + '<span class="bolded-red">' + seq[11] + '</span>' + seq[6:11]
        sorted_val = sorted(tbl[row_key],reverse=True,key=lambda x:abs(x[1]))
        for row_val in sorted_val:
            if filteropt == 1:
                pbmname = row_val[2]
            else:
                pbmname = row_val[3]
            if pbmname == 'None':
                ingenes_str = ""
            else:
                ingenes_str = "\"" + ",".join([gene for gene in pbmtohugo[pbmname] if gene in gene_names]) + "\""
            if filteropt == 1:
                datavalues.append([row_key[0],seq[0:11],(seq[0:5] + seq[11] + seq[6:11]),row_val[0],row_val[1],pbmname,ingenes_str])
                #csv_ret+=("{},{},{},{},{:.3f},{},{}\n".format(row_key[0],seq[0:11],(seq[0:5] + seq[11] + seq[6:11]),row_val[0],row_val[1],pbmname,ingenes_str))
            else:
                datavalues.append([row_key[0],seq[0:11],(seq[0:5] + seq[11] + seq[6:11]),row_val[0],Decimal(row_val[1]),row_val[2],pbmname,ingenes_str])
                #csv_ret+=("{},{},{},{},{:.3e},{},{},{}\n".format(row_key[0],seq[0:11],(seq[0:5] + seq[11] + seq[6:11]),row_val[0],Decimal(row_val[1]),row_val[2],pbmname,ingenes_str))
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
                        if abs(row_val[1]) > abs(maintbl[row_key][least_idx][1]):
                            del maintbl[row_key][least_idx]
                            maintbl[row_key].append(row_val)
            else: # filteropt == 2 -- p-value
                for row_key in ddict:
                    maintbl[row_key].extend(ddict[row_key])
    return format2tbl(maintbl,gene_names,filteropt)

@celery.task()
def delete_file(filename):
    '''
    this simple function is used to delete user file after USER_DATA_EXPIRY
    seconds
    '''
    if os.path.exists(filename):
        os.remove(filename)
        print("Deleted: %s"%filename)
    else:
        print("%s doesn't exist for deletion"%filename)

#==========================================================

'''
intbl: preprocessed table
filteropt: 1 for highest t-val, 2 for p-val cutoff
filterval: # TFs for opt 1 and p-val cutoff for opt 2
'''
@celery.task(bind=True)
def do_prediction(self,intbl,selections,gene_names,filteropt=1,filterval=1):
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

    sharedlist = mp.Manager().list()
    async_pools = [pool.apply_async(predict, (preds[i],intbl,sharedlist,filteropt,filterval)) for i in range(0,len(preds))]

    # run the job, update progress bar
    total = len(predfiles)
    while not all([p.ready() for p in async_pools]):
        ready_sum = len(sharedlist)
        time.sleep(2) # super important to avoid checking every loop
        self.update_state(state='PROGRESS',
                          meta={'current': ready_sum, 'total': total, 'status': 'Processing input data...'})

    res = [p.get() for p in async_pools]
    self.update_state(state='PROGRESS',
                          meta={'current': ready_sum, 'total': total, 'status': 'post-processing'})
    print("Terminate all children process..")
    pool.terminate() # terminate to kill all child processes !!! Like.. super important,
                     # to avoid memory leak, seriously...
    colnames,datavalues = postprocess(res,gene_names,filteropt,filterval)

    if db.exists(self.request.id):
        reqdict = db.hgetall(self.request.id)
    else: # not likely to happen though
        reqdict = {}
    reqdict['rescol'] = colnames
    reqdict['resval'] = datavalues
    if db.exists(self.request.id):
        db.delete(self.request.id)
    db.hmset(self.request.id,reqdict)
    db.expire(self.request.id, app.config['USER_DATA_EXPIRY'])

    csv_path = "%s.csv"%self.request.id # delete this
    #with open(app.config['UPLOAD_FOLDER'] + csv_path,'w') as f:
    #    f.write(csv_ret)

    # https://stackoverflow.com/questions/24577349/flask-download-a-file
    return {'current': len(sharedlist), 'total': len(predfiles), 'status': 'Task completed!',
            'result': 'done', 'taskid': self.request.id,
            'time':(time.time()-start_time)} # -- somehow cannot do jsonify(postproc)
