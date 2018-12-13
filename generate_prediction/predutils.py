import argparse
import pandas as pd
import os

def get_gapdata(pbmpath):
    '''
    print gappos and then gapsize for shell script to parse
    '''
    pbmname = os.path.splitext(os.path.basename(pbmpath))[0]
    db = "resource/upbm_gap_params.csv"
    df = pd.read_csv(db)
    row = df.loc[df['upbm_filenames'] == pbmname]
    print(row['gappos'].item())
    print(row['gapsize'].item())

#python3 predutils.py -p resource/data_upbm_selected_good_quality_final.csv
def gapparams_from_csv(csvpath):
    df = pd.read_csv(csvpath)
    datacsv = []
    processed = []
    for idx,row in df.iterrows():
        upbmname = os.path.splitext(row['upbm_filenames'])[0]
        if row['upbm_filenames'] in processed:
            continue
        else:
            processed.append(upbmname)
        rowdict = {}
        rowdict['upbm_filenames'] = upbmname
        rowdict['gapmodel'] = row['best']
        if row['best'] == "ungapped":
            rowdict['gapsize'] = 0
            rowdict['gappos'] = 0
        elif row['best'].startswith("gap"):
            gaps = row['best'][len("gap"):].split("p")
            rowdict['gapsize'] = gaps[0]
            rowdict['gappos'] = gaps[1]
        datacsv.append(rowdict)
    pd.DataFrame(datacsv).to_csv("upbm_gap_params.csv",index=False,columns=['upbm_filenames','gappos','gapsize','gapmodel'])

#python3 predutils.py -g "test-in/Mus_musculus|NA|Unpublished|Zfp24.txt"
if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Generate predicitons for all possible k-mer permutations.')
    parser = argparse.ArgumentParser(description='Utilities to help generate prediction.')
    parser.add_argument('-p','--parse', type=str, help="Parse an input csv to generate parameters for the gap models")
    parser.add_argument('-g','--gappbm', type=str, help="Get the gap parameters for a pbmname")
    args = parser.parse_args()

    if args.parse is not None:
        gapparams_from_csv(args.parse)
    elif args.gappbm is not None:
        get_gapdata(args.gappbm)
