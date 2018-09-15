# app/__init__.py

from flask import Flask
from celery import Celery
import redis

# Initialize the app
app = Flask(__name__, instance_relative_config=True)
app.secret_key = "super secret key"

# Load the config file
app.config.from_object('config')

celery = Celery(app.name, backend=app.config['CELERY_RESULT_BACKEND'], broker=app.config['CELERY_BROKER_URL'])
celery.conf.update(app.config)

db = redis.Redis('localhost', decode_responses=True) # decode->make it in utf8,TODO

# Load the views
from app import views
from app import test
