# run.py

from app import app,db
#import app.model.redispubsub as redispubsub

if __name__ == '__main__':
    #pubsub = redispubsub.run_pubsub_thread(db)
    app.run(use_reloader=False,debug=True) # (debug=True,host='0.0.0.0')

    #use_reloader: somehow make this file to be called twice in debug mode
