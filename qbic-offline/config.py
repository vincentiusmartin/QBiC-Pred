import configparser
import importlib.util
import pandas as pd

import os

def get_family_map(topbm_mapping_path):
    with open(topbm_mapping_path,'r') as f:
        family_map = {}
        for line in f:
            key,val = line.strip().split("->")
            valmap = {z[0]:z[1] for z in (y.split(":") for y in (x for x in val.split(";")))} # generator
            family_map[key] = valmap
    return family_map

def parse_hugo_name_mapping(filepath):
    df = pd.read_csv(filepath,sep=" ")
    mapping = pd.Series(df.hugo_id.values,index=df.hugo_name).to_dict()
    return mapping

def dictfamily2genedict(dictlist):
    result = {}
    for d in dictlist:
        new_d = {k:d[k].split(",") for k in d}
        result.update(new_d)
    return result

config = configparser.ConfigParser()
if os.path.exists("config.ini"):
    config.read('config.ini')
else:
    raise FileNotFoundError
    sys.exit(1)

if config["General Conf"]["PCOUNT"] == "cpu.count":
    PCOUNT = os.cpu_count()
else:
    PCOUNT = int(config["General Conf"]["PCOUNT"])

''' [Directory Setting] '''
PREDDIR = config["Directory Setting"]["PREDDIR"]
CHRDIR = config["Directory Setting"]["CHRDIR"]
ESCORE_DIR = config["Directory Setting"]["ESCORE_DIR"]
PBM_HUGO_MAPPING = config["Directory Setting"]["PBM_HUGO_MAPPING"]
DBD_PBM_MAPPING = get_family_map(config["Directory Setting"]["HUGO_PBM_MAPPING"])
HUGO_PBM_MAPPING = dictfamily2genedict([DBD_PBM_MAPPING[key] for key in DBD_PBM_MAPPING])
