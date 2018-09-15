# config.py

# Enable Flask's debugging features. Should be False in production
DEBUG = True

# Celery related configurations
CELERY_BROKER_URL = 'redis://localhost:6379/0'
CELERY_RESULT_BACKEND = 'redis://localhost:6379/0' # important to track result
#CELERYD_TASK_SOFT_TIME_LIMIT = 60

# Flask configs
MAX_CONTENT_LENGTH = 32 * 1024 * 1024 # max upload size

PCOUNT = 1 #os.cpu_count()
PREDDIR = "/Users/vincentiusmartin/Research/MutationPredictor/tfbc-website/preddir" #"/usr/project/xtmp/vmartin/pred"
UPLOAD_FOLDER = '/tmp/'
CHRDIR = "/Users/vincentiusmartin/Research/MutationPredictor/tfbc-website/chromosomes"
MAPPING_FILE = "resource/tf-mapping.txt"
