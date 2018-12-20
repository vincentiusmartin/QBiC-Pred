import sys
sys.path.insert(0, '..')

from flask import jsonify,request,url_for

from app import app,db

import json

@app.route('/recent', methods=['GET'])
def get_recent_jobs():
    recents = []
    for key in request.cookies.keys():
        if key.startswith("qbic_recents:"):
            id = key.split(":",1)[1]
            if db.exists(id):
                recents.append(key)
    if recents:
        rj_urls = [[request.cookies.get(job_key),url_for('process_request',job_id=job_key.split(":",1)[1])] for job_key in recents]
        return json.dumps(rj_urls)
    else:
        return json.dumps({})
