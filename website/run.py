# run.py

from app import app

if __name__ == '__main__':
    app.run(host='0.0.0.0',use_reloader=False,debug=True) # (debug=True,host='0.0.0.0')

    #use_reloader: somehow make this file to be called twice in debug mode
