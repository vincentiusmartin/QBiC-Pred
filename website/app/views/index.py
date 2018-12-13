from flask import render_template,request,jsonify,make_response,url_for
from werkzeug.utils import secure_filename
from app import app,db
import os
import uuid

import pandas as pd

import app.controller.celerytask as celerytask

@app.route('/', methods=['GET', 'POST'])
def index():
    #session.permanent = True
    #session.clear() -- need to limit the amount of session somewhere
    return render_template("index.html")

@app.route('/submitpredfile', methods=['POST'])
def submit_pred_upload():
    if 'predupload-file' not in request.files:
        return jsonify({'Message':'no input file part'}), 500
    file = request.files['predupload-file']
    filename = secure_filename(file.filename)
    filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
    file.save(filepath)
    df = pd.read_csv(filepath,dtype=str)
    os.remove(filepath)

    rand_id = str(uuid.uuid4())

    cols = list(df.columns.values)
    if "z_score" in cols:
        filteropt = 1
    else:
        filteropt = 2
    filterval = "-"
    genes_selected = df.TF_gene.unique()
    datavalues = df.to_dict('records')

    celerytask.savetoredis(rand_id,cols,datavalues,app.config['UPLOAD_PRED_EXPIRY'])

    session_info = {"parent_id":"uploadpred",
                    "task_id":rand_id,
                    "filename":"",
                    "genes_selected":genes_selected,
                    "filteropt":filteropt,
                    "filterval":filterval,
                    "chrver":""}
    if db.exists(rand_id):
        db.delete(rand_id)
    db.hmset(rand_id,session_info)
    db.expire(rand_id, app.config['UPLOAD_PRED_EXPIRY'])

    resp = make_response(jsonify({}), 202, {'Location': url_for('process_request',job_id=rand_id)})
    return resp
