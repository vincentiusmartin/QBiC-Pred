import configparser
import importlib.util

def import_from_file(filepath):
    # https://stackoverflow.com/questions/67631/how-to-import-a-module-given-the-full-path
    spec = importlib.util.spec_from_file_location("", filepath)
    foo = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(foo)
    return foo

config = configparser.ConfigParser()
config.read('qbic-conf.ini')

''' [Debug] '''
# Enable Flask's debugging features. Should be False in production
DEBUG = config["Debug"]["DEBUG"]

''' [Celery Conf] '''
# Celery related configurations
CELERY_BROKER_URL = config["Celery Conf"]["CELERY_BROKER_URL"]
CELERY_RESULT_BACKEND = config["Celery Conf"]["CELERY_RESULT_BACKEND"]

''' [Flask Conf] '''
MAX_FILE_LENGTH =  int(config["Flask Conf"]["MAX_FILE_LENGTH"]) # CONTENT

if config["Flask Conf"]["PCOUNT"] == "cpu.count":
    PCOUNT = os.cpu_count()
else:
    PCOUNT = int(config["Flask Conf"]["PCOUNT"])

''' [Directory Setting] '''
PREDDIR = config["Directory Setting"]["PREDDIR"]
UPLOAD_FOLDER = config["Directory Setting"]["UPLOAD_FOLDER"]
CHRDIR = config["Directory Setting"]["CHRDIR"]
ESCORE_DIR = config["Directory Setting"]["ESCORE_DIR"]
PBM_HUGO_MAPPING = config["Directory Setting"]["PBM_HUGO_MAPPING"]
HUGO_PBM_MAPPING = config["Directory Setting"]["HUGO_PBM_MAPPING"]
GAP_FILE = config["Directory Setting"]["GAP_FILE"]
STATIC_EXAMPLE_DIR = config["Directory Setting"]["STATIC_EXAMPLE_DIR"]
INPUT_EXAMPLE_DICT = import_from_file(config["Directory Setting"]["INPUT_EXAMPLE_LIST"]).examples


''' [User Session] '''
USER_DATA_EXPIRY = int(config["User Session"]["USER_DATA_EXPIRY"])
UPLOAD_PRED_EXPIRY = int(config["User Session"]["UPLOAD_PRED_EXPIRY"])
