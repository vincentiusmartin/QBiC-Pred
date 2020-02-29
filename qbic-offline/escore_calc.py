import utils
import argparse

import pandas as pd
from timeit import default_timer as timer

def isbound_escore_18mer_list(seq18mer_list,escore_path,short2long_map="",spec_ecutoff=0.35,nonspec_ecutoff=0.4):
    #eshort_path = "%s/%s_escore.txt" % (escore_dir,pbm_name)
    # TODO: avoid IO, maybe using global var?
    #short2long_map = "%s/index_short_to_long.csv" % (escore_dir)

    #  making escore using the concise version, use integer indexing
    if short2long_map:
        with open(escore_path) as f:
            eshort = [float(line) for line in f]

        with open(short2long_map) as f:
            next(f)
            emap = [int(line.split(",")[1])-1 for line in f]

        elong = [eshort[idx] for idx in emap]
    else: # making escore using normal escore file
        df = pd.read_csv(escore_path, sep="\t")
        d1 = pd.Series(df["E-score"].to_numpy(),index=df["8-mer"]).to_dict()
        d2 = pd.Series(df["E-score"].to_numpy(),index=df["8-mer.1"]).to_dict()
        elong = {**d1, **d2}

    bound_list = []
    for seq in seq18mer_list:
        wild = seq[:-1]
        mut = seq[:8] + seq[-1] + seq[9:-1]
        if short2long_map:
            change = "%s>%s" % (utils.isbound_escore(wild,elong,bsite_cutoff=spec_ecutoff,nbsite_cutoff=nonspec_ecutoff),
                              utils.isbound_escore(mut,elong,bsite_cutoff=spec_ecutoff,nbsite_cutoff=nonspec_ecutoff))
        else:
            change = "%s>%s" % (utils.isbound_escore_8merdict(wild,elong,bsite_cutoff=spec_ecutoff,nbsite_cutoff=nonspec_ecutoff),
                              utils.isbound_escore_8merdict(mut,elong,bsite_cutoff=spec_ecutoff,nbsite_cutoff=nonspec_ecutoff))
        bound_list.append({"wild":wild, "mutant":mut, "change":change})

    return bound_list

if __name__ == "__main__":
    """
    inputfile = "testing_resources/single_mutations_sample.tsv"
    eshort_path = "testing_resources/E2F1_8mers_escore_compact.txt"
    efull_path = "testing_resources/E2F1_8mers_escore_full.txt"
    eshort2long_path = "/Users/vincentiusmartin/Research/MutationPredictor/QBiC-Pred/escore/escore/index_short_to_long.csv"
    chrver = "hg19"
    filetype = "mut"
    """
    import qbic # not a necessary import for isbound, so save an import conditionally

    parser = argparse.ArgumentParser(description = 'TF Mutation Predictions')
    parser.add_argument('-i', '--inputfile', action="store", dest="inputfile", type=str,
                        help='Input mutation file in .vcf, .tsv, .csv, or .txt format.')
    parser.add_argument('-e', '--escorepath', action="store", dest="escorepath", type=str,
                        help="""
                        Input e-score file, by default the full e-score file is expected. Othrwise,
                        if the concise version is used then indexmap needs to be set.
                        """)
    parser.add_argument('-m', '--indexmap', action="store", dest="indexmap", type=str, default="",
                        help="""Escore concise to long map. If this option is set then the program will
                        assume that the concise Escore file is used.""")
    parser.add_argument('-c', '--chrver', action="store", dest="chrver", type=str, default="hg19",
                        help="Chromosome version to use.")
    parser.add_argument('-t', '--ftype', action="store", dest="ftype", type=str, default="",
                        help="""File type can specify: vcf, icgc, customseq, or mut.
                        If it is not specified then file extension will be checked.""")
    args = parser.parse_args()


    # concise escore: python3 escore_calc.py -i testing_resources/single_mutations_sample.tsv -e testing_resources/E2F1_8mers_escore_compact.txt -m /Users/vincentiusmartin/Research/MutationPredictor/QBiC-Pred/escore/escore/index_short_to_long.csv -t mut
    # full escore: python3 escore_calc.py -i testing_resources/single_mutations_sample.tsv -e testing_resources/E2F1_8mers_escore_full.txt -t mut

    allseq = [arr[2] for arr in qbic.inittbl(args.inputfile, args.chrver, filetype = args.ftype)]
    test_start = timer()
    escore_res = isbound_escore_18mer_list(allseq, args.escorepath, args.indexmap) # eshort2long_path
    df = pd.DataFrame(escore_res)
    df.to_csv("escore_result.csv")
    test_end = timer()
    total_time = test_end-test_start
    print("Total E-score calculation time {:.2f}secs".format(total_time))
