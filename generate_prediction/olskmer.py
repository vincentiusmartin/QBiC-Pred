import pandas as pd
import statsmodels.formula.api as sm
import numpy as np
import itertools
import time
import os,sys
import argparse

#sys.path.append('../localpackage')
import bio

from sklearn import linear_model

nucleotides = ['A','C','G','T']

def adjustscr(score,shift=1000):
    minscr = score[score.idxmin()]
    score_frame = pd.DataFrame(np.log(score - minscr + shift),columns=['score'])
    return score_frame

def read_pbm(filename,kmer):
    tbl = pd.read_csv(filename,names=['score','sequence'],delim_whitespace=True) #score,sequence ..
    score = adjustscr(tbl['score'],1000)
    seqbin = [bio.seqtoi(x) for x in tbl['sequence']]
    oligfreq = bio.nonr_olig_freq(seqbin,kmer) #tbl['sequence'].tolist()
    return pd.concat([score,oligfreq],axis=1)

def print_full(x):
    y = pd.DataFrame(x)
    #y = pd.DataFrame([bio.itoseq(z) for z in x])
    pd.set_option('display.max_rows', len(y))
    y.columns = [bio.itoseq(x) for x in y.columns]
    print(y)
    pd.reset_option('display.max_rows')

if __name__ == "__main__": # python olskmer.py $filein $output -k $kmer -d $chunk
    parser = argparse.ArgumentParser(description='Generate predicitons for all possible k-mer permutations.')
    parser.add_argument('input', type=str, nargs=2, help="Input for the script, format: <tf-path> <out-path>")
    parser.add_argument('-k','--kmer', type=int, help="The value of k that we want", default=6)
    parser.add_argument('-d','--div', type=int, help="Sometimes data can be too big, div is the number of chunks required", default=32)
    args = parser.parse_args()

    tfpath = args.input[0]
    outpath = args.input[1]

    filename = os.path.splitext(os.path.basename(tfpath))[0]
    print("Got input file: " + filename)

    bio.gen_noreversed_kmer(args.kmer)  
    start_time = time.time()
    df = read_pbm(tfpath,args.kmer)  
    print("---Time to read pbm file: %s seconds ---" % (time.time() - start_time)) #5minutes
    lm = sm.OLS(df['score'],df.drop('score',axis=1)).fit()
    print("--- Finish training models: %s seconds ---" % (time.time() - start_time)) 

    #mutated context
    mutated = dict() 
    for base in bio.nucleotides: # each is sequence of length 11 with the mutated in the middle
        mutated[base] = [bio.insert_pos(x,base,args.kmer-1) for x in bio.seq_permutation(2*(args.kmer-1))]  
    # we have 1048576 combinations from 4**10
    
    chunk = len(bio.nucleotides)**((args.kmer-1)*2) // args.div
    #print("Chunk size: " + str(chunk))
    output_all = pd.DataFrame(columns=['dna_seq','diff','t'])
    for i in range(0,args.div):
        print("Processing chunk-"+str(i))
        mutated_part = dict()

        count = dict() 
        # count oligonucleotides in all k-mer combinations of sequence length k with mutated in the middle.
        # row = all sequences input, col = all possible permutations
        # Therefore we have 4**2k-1 rows (i.e. combinations) 
        for base in bio.nucleotides:          
            mutated_part[base] = mutated[base][i*chunk:(i+1)*chunk] #16384
            count[base] = bio.nonr_olig_freq(mutated_part[base],args.kmer) #count the frequency of kmer in the mutation ##save this?
            # for k == 4: 4096 x 136
        #dff = pd.DataFrame(count['A'])
        #dff['xyz'] = [bio.itoseq(ddd) for ddd in mutated_part['A']]        
        #print_full(dff) # TTTATTT = 4095

        diff_count = dict()
        diff = dict()
        for b1 in nucleotides:
            single = []
            for b2 in nucleotides:
                if b1 != b2:
                    single.append(count[b2] - count[b1])
            diff_count[b1] = pd.concat(single,axis=0,ignore_index=True)
            diff[b1] = np.dot(diff_count[b1],lm.params) #dim: (49152, 2080) (2080,) .. diff_count == c' in the paper?
        #print("{} {}".format(diff_count[b1].shape,lm.params.shape))        
        del count        
        
        sd_diff = dict()
        t = dict() #for t-test
        for base in nucleotides: #math.sqrt((diff_count[base].transpose() * np.dot(lm.cov_params(),diff_count[base].transpose())).sum(axis=0))
            sd_diff[base] = (diff_count[base].transpose() * np.dot(lm.cov_params(),diff_count[base].transpose())).sum(axis=0).apply(np.sqrt)                     
            t[base] = diff[base] / sd_diff[base]
        
        dna_seq = []
        diff_all = []
        t_all = []
        for b1 in nucleotides:
            diff_all += diff[b1].tolist()
            t_all += t[b1].tolist()
            for b2 in nucleotides:
                if b1 != b2:
                    dna_seq += [ ((x << 2) | bio.nucleotides[b2]) for x in mutated_part[b1]]
                    
        newout = pd.DataFrame({'dna_seq':[bio.itoseq(x) for x in dna_seq],
                                        'diff':diff_all,
                                        't':t_all},columns=['dna_seq','diff','t'])
        output_all = output_all.append(newout,ignore_index = True)

    na_entries = pd.DataFrame(columns=['dna_seq','diff','t'])   
    for base in nucleotides:
        na_entries = na_entries.append(pd.DataFrame({
            'dna_seq': [bio.itoseq((bio.insert_pos(x,base,args.kmer-1) << 2) | bio.nucleotides[base]) for x in bio.seq_permutation(2*(args.kmer-1))],
            'diff':np.nan,
            't':np.nan
        },columns=['dna_seq','diff','t']))
    output_all = output_all.append(na_entries,ignore_index = True).sort_values(['dna_seq'],ascending=True) #replace(np.nan, 'NaN', regex=True)
    output_all.to_csv("{}/prediction{}mer.{}.csv".format(outpath,args.kmer,filename),columns=['diff','t'],sep=' ',index=None,float_format="%.5f")
    
    print("--- Total time: %s seconds ---" % (time.time() - start_time))
#4372.6606secs

