'''
Output type: 1 for largest absolute z-score and 2 for p-val threshold
'''

examples = \
{
    'icgc-breast-cancer-mutations':{
        'inputfile':'icgc-example.tsv',
        'tfs':['MYC','E2F','FOXA1'],
        'genomever':"hg19",
        'outputtype':1
    },
    'MAFK-ChIP-seq-allele-specific-binding-variants':{
        'inputfile':'mafk-asb-example.tsv',
        'tfs':['MAFK'],
        'genomever':"hg19",
        'outputtype':2
    },
}
