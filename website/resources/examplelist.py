'''
Output type: 1 for largest absolute z-score and 2 for p-val threshold
'''

examples = \
{
    'icgc-breast-cancer-mutations':{
        'inputfile':'icgc-example.tsv',
        'tfs':['E2F1','FOXA1','MYC'],
        'genomever':"hg19",
        'outputtype':1
    },
    'MAFK-ChIP-seq-allele-specific-binding-variants':{
        'inputfile':'mafk-asb-example.vcf',
        'tfs':['MAFK'],
        'genomever':"hg19",
        'outputtype':2
    },
}
