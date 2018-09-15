
from celery import Celery,chain
from app import app,celery

import time,sys,gzip
import pandas as pd
import billiard as mp #multiprocessing substitute to enable daemon
import scipy.stats

import os

def get_chrom(cfile):
    with gzip.open(cfile,'rb') as f:
        next(f)
        chrom = f.read().decode('utf-8').replace('\n','')
    return chrom

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

# https://stackoverflow.com/questions/2130016/splitting-a-list-into-n-parts-of-approximately-equal-length
def chunkify(lst,n):
    return [lst[i::n] for i in range(n)]

# =====================================================

# input: table
@celery.task(bind=True)
def inittbl(self,filename,cpath):
    self.update_state(state='PROGRESS',
                  meta={'current': 0, 'total': 1, 'status': 'Preprocessing input...'})
    kmer = 6    # TODO: handle separator csv/tsv
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
    for cidx in [str(a) for a in range(1,23)] + ['X','Y']:
        self.update_state(state='PROGRESS',
                  meta={'current': 0, 'total': 1, 'status': 'Preprocessing input for chromosome {}...'.format(cidx)})
        if cidx not in dataset:
            continue
        print("Iterating dataset for chromosome {}...".format(cidx))
        chromosome = get_chrom(cpath + "/chr." + str(cidx) + '.fa.gz')
        for idx,row in dataset[cidx].iterrows():
            pos = row['chromosome_start'] - 1
            if row['mutated_from_allele'] != chromosome[pos]:
                sys.exit("Found mismatch in the mutation: \n{}".format(row))
            seq = chromosome[pos-kmer+1:pos+kmer] + row['mutated_to_allele'] #-5,+6
            result.append([idx,seq,seqtoi(seq),0,0,"-"]) #rowidx,seq,val,diff,t,pbmname

    result = sorted(result,key=lambda result:result[0])
    print("Time to preprocess: {:.2f}secs".format(time.time()-start))
    return result

#==================================== Prediction Part ====================================

# TODO: remove sharedlist if not needed anymore
'''
for the container list, key is a tuple of: (rowidx,sequence,seqidx)
and each element in value is a list if: [diff,z-score,pbmname]
'''
def predict(predlist,dataset,sharedlist,filteropt=1,filterval=1):
    buggedtf = 0
    if filteropt == 1:
        container = {tuple(row):[[0,0,"None"]] for row in dataset} # diff,t,pbmname -- 0,0,"-"
    else:
        container = {tuple(row):[] for row in dataset}

    # iterate for each transcription factor
    for i in range(0,len(predlist)):
        start = time.time()
        tfname = '.'.join(map(str,predlist[i].split(".")[1:-1]))
        print("Processing ",tfname)
        with open(predlist[i], 'r') as f:
            tflist = pd.read_csv(f, delimiter=' ').round(5).values.tolist()
        if len(tflist) < 4**12:
            print("Skip %s since it has less rows than 4**12" % tfname)
            buggedtf += 1
            continue
        for row_key in container:
            seqidx = row_key[2]
            if filteropt == 1:
                if len(container[row_key]) >= filterval:
                    least_idx = min(enumerate(container[row_key]),key=lambda x:abs(x[1][1]))[0]
                    if abs(tflist[seqidx][1]) > abs(container[row_key][least_idx][1]):
                        del container[row_key][least_idx]
                        container[row_key].append(tflist[seqidx] + [tfname])
                else:
                    container[row_key].append(tflist[seqidx] + [tfname])
            else: # filteropt == 2
                zscore = tflist[seqidx][1]
                pval = scipy.stats.norm.sf(abs(zscore))*2
                if pval <= filterval:
                    #print(list(row_key)+tflist[seqidx][:-1] + [pval,tfname])
                    container[row_key].append(tflist[seqidx][:-1] + [pval,tfname])
        sharedlist.append(tfname) # TODO: might not be needed
        print("Total running time for {}: {:.2f}secs".format(tfname,time.time()-start))

    '''if filteropt == 2:
        for key in container:
            if not container[key]:
                container[key] = [0,0,"None"]'''

    # remove seqidx as it is not needed anymore
    newcontainer = {}
    for row_key in container:
        newcontainer[row_key[:-1]] = container[row_key]

    return newcontainer

def format2tbl(tbl,filteropt=1):
    if filteropt == 1:
        csv_ret = "rowidx,wild,mutant,diff,z-score,pbmname\n"
        metrics = 'z-score'
    else: #filteropt == 1:
        csv_ret = "rowidx,wild,mutant,diff,p-value,pbmname\n"
        metrics = 'p-value'
    sorted_key = sorted(tbl.keys())
    for row_key in sorted_key:
        if not tbl[row_key]: # probably empty row
            continue
        seq = row_key[1]
        wild = seq[0:5] + '<span class="bolded-red">' + seq[5] + '</span>' + seq[6:11]
        mut = seq[0:5] + '<span class="bolded-red">' + seq[11] + '</span>' + seq[6:11]
        sorted_val = sorted(tbl[row_key],reverse=True,key=lambda x:abs(x[1]))
        for row_val in sorted_val:
            csv_ret+=("{},{},{},{},{:.3f},{}\n".format(row_key[0],seq[0:11],(seq[0:5] + seq[11] + seq[6:11]),row_val[0],row_val[1],row_val[2]))
    '''
    tblret.append({'rowidx':row_key[0],'wild-type':wild,'mutant':mut,'diff':row_val[0],metrics:"%.5f"%row_val[1],'pbmname':row_val[2]})
    #tblret += '{\'rowidx\':%s,\'wild-type\':%s,\'mutant\':%s,\'diff\':%s,\'%s\':%.5f,\'pbmname\':%s},' % (row_key[0],wild,mut,row_val[0],metrics,row_val[1],row_val[2])
    tblret.append({'rowidx':row_key[0],'wild-type':wild,'mutant':mut,'diff':row_val[0],metrics:"%.5f"%row_val[1],'pbmname':row_val[2]})
    tblret = tblret[:-1]
    '''
    return csv_ret

def postprocess(datalist,filteropt=1,filterval=1):
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
    return format2tbl(maintbl,filteropt)

#==========================================================

'''
intbl: preprocessed table
filteropt: 1 for highest t-val, 2 for p-val cutoff
filterval: # TFs for opt 1 and p-val cutoff for opt 2
'''
@celery.task(bind=True)
def do_prediction(self,intbl,selections,filteropt=1,filterval=1):
    start_time = time.time()

    #while not inittask.ready():
    #    time.sleep(1)
    #intbl = inittask.get()

    # move the comment here for testing
    pool = mp.Pool(processes=app.config['PCOUNT'])
    predfiles = [app.config['PREDDIR'] + "/" + s for s in selections] # os.listdir(preddir)
    preds = chunkify(predfiles,app.config['PCOUNT']) # chunks the predfiles for each process

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
    csv_ret = postprocess(res,filteropt,filterval)

    csv_path = "%s.csv"%self.request.id
    with open(app.config['UPLOAD_FOLDER'] + csv_path,'w') as f:
        f.write(csv_ret)

    # https://stackoverflow.com/questions/24577349/flask-download-a-file
    return {'current': len(sharedlist), 'total': len(predfiles), 'status': 'Task completed!',
            'result': 'done','csvlink':csv_path,
            'time':(time.time()-start_time)} # -- somehow cannot do jsonify(postproc)
    '''
    i = 0
    while i < len(intbl):
        print("lalala " + str(intbl[i]))
        time.sleep(1)
        self.update_state(state='PROGRESS',
                          meta={'current': i, 'total': len(intbl), 'status': 'progressing'})
        i+= 1

    return {'current': i, 'total': len(intbl), 'status': 'Task completed!',
            'result': [{0:'rowidx',1:'sequence'},{'rowidx':200,'sequence':'TATGCG'}]} # dummy
    '''
