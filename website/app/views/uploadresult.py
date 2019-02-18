from flask import render_template,request,jsonify,make_response,url_for
from werkzeug.utils import secure_filename
from app import app,db
import os
import uuid

import pandas as pd
from pandas.core.groupby.groupby import DataError

import app.controller.celerytask as celerytask
import app.controller.utils as utils

@app.route('/uploadresult', methods=['GET', 'POST'])
def upload_result():
    #session.permanent = True
    #session.clear() -- need to limit the amount of session somewhere
    return render_template("uploadresult.html")

def prepare_predfile(request):
    if 'predupload-file' not in request.files:
        return 'error', 'no input file part'
    file = request.files['predupload-file']
    filename = secure_filename(file.filename)
    filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
    file.save(filepath)

    try:
        df = pd.read_csv(filepath,dtype=str)
    except:
        return 'error', 'input is not supported'
    utils.delete_file(filepath) #delete once read
    check_cols = set(["row","wild","mutant","diff","z_score","p_value","TF_gene","binding_status","gapmodel","pbmname"])
    df_cols = set(df.columns)
    if not check_cols.issubset(df_cols):
        return 'error', 'could not find all required fields'
    return 'success',df

@app.route('/submitpredfile', methods=['POST'])
def submit_pred_upload():
    status,message = prepare_predfile(request)
    if status == "error":
        return jsonify({'Message':message}), 500

    df = pd.DataFrame(message) # if success then message is the dataframe
    rand_id = str(uuid.uuid4())

    cols = list(df.columns.values)
    if "z_score" in cols:
        filteropt = 1
    else:
        filteropt = 2
    filterval = "-"
    genes_str = ",".join(list(df.TF_gene))
    genes_selected = list(set(genes_str.split(",")))
    datavalues = df.to_dict('records')

    celerytask.savetoredis(rand_id,cols,datavalues,app.config['UPLOAD_PRED_EXPIRY'])

    session_info = {"parent_id":"uploadpred",
                    "task_id":rand_id,
                    "filename":"-",
                    "genes_selected":genes_selected,
                    "filteropt":filteropt,
                    "filterval":filterval,
                    "chrver":"-",
                    "spec_escore_thres":"-",
                    "nonspec_escore_thres":"-"}
    if db.exists(rand_id):
        db.delete(rand_id)
    db.hmset(rand_id,session_info)
    db.expire(rand_id, app.config['UPLOAD_PRED_EXPIRY'])

    resp = make_response(jsonify({}), 202, {'Location': url_for('process_request',job_id=rand_id)})
    return resp
