'''
Output type: 1 for largest absolute z-score and 2 for p-val threshold
'''

examples = \
{
    'example1':{
        'inputfile':'icgc-example.tsv',
        'tfs':['TfX','E2F'],
        'genomever':"hg38",
        'outputtype':1
    },
    'example2':{
        'inputfile':'icgc-example.tsv',
        'tfs':['TfX'],
        'genomever':"hg19",
        'outputtype':2
    },
}
