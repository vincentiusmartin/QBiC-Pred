import app.controller.utils as utils

# this assume that inputs are already correct
with open("app/static/files/mutated-seqs.txt") as f:
    result = []
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
    print(result)
