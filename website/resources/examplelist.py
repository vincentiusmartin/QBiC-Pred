'''
Output type: 1 for largest absolute z-score and 2 for p-val threshold
'''

examples = \
{
    'MAFK ChIP-seq Allele Specific Binding Variants':{
        'inputfile':'mafk-asb-example.vcf',
        'tfs':['MAFK'],
        'genomever':"hg19",
        'outputtype':2
    },
    'ICGC Breast Cancer Mutations - Small':{
        'inputfile':'icgc-example-short.tsv',
        'tfs':['E2F1','FOXA1','MYC'],
        'genomever':"hg19",
        'outputtype':1
    },
    'ICGC Breast Cancer Mutations - Large':{
        'inputfile':'icgc-example.tsv',
        'tfs':['E2F1','FOXA1','MYC'],
        'genomever':"hg19",
        'outputtype':1
    }
}
