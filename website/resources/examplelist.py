'''
Output type: 1 for largest absolute z-score and 2 for p-val threshold
'''

examples = \
{
    'MAFK ChIP-seq Allele Specific Binding Variants':{
        'inputfile':'QBiC-vcf-example-MAFK-ASB-variants.vcf',
        'tfs':['MAFK'],
        'genomever':"hg19",
        'outputtype':2
    },
    'ICGC Breast Cancer Mutations - Small':{
        'inputfile':'QBiC-icgc-example-breast-cancer-mutations-small.tsv',
        'tfs':['E2F1','FOXA1','MYC'],
        'genomever':"hg19",
        'outputtype':1
    },
    'ICGC Breast Cancer Mutations - Large':{
        'inputfile':'QBiC-icgc-example-breast-cancer-mutations-large.tsv',
        'tfs':['E2F1','FOXA1','MYC'],
        'genomever':"hg19",
        'outputtype':1
    }
}
