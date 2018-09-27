# views.py

import os

from flask import Flask,flash,request,redirect,url_for,render_template,send_from_directory,jsonify,session,make_response
from celery import Celery,chain
from celery.task.control import revoke

from app import app,celery,db

from werkzeug.utils import secure_filename
from decimal import Decimal

import json,ast
import pandas as pd

import config

import app.celerytask as celerytask

ALLOWED_EXTENSIONS = set(['csv','tsv'])
#app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

# view specific utils
def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

#===============JOB STATUS/ RESULT PAGE RELATED=================

@app.route('/files/<path:filename>')
def get_file(filename):
    """Download a file."""
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename,as_attachment=True)

def filter_search(infields,keyword,type):
    if type == "All":
        if any(keyword in str(field) for field in infields):
            return True
        else:
            return False
    return True

@app.route('/getrescol/<task_id>',methods=['GET'])
def get_res_col(task_id):
    with open("%s%s.csv"%(app.config['UPLOAD_FOLDER'],task_id)) as f:
        cols = [{"title":title} for title in f.readline().strip().split(",")]
    return jsonify(cols)

@app.route('/getrestbl/<task_id>',methods=['GET'])
def get_res_tbl(task_id):
    # {\"0\":\"rowidx\",\"1\":\"wild-type\",\"2\":\"mutant\",\"3\":\"diff\",\"4\":\"z-score\",\"5\":\"pbmname\"},
    # Get all mandatory informations we need
    draw = int(request.args['draw']) # not secure # TODO: make it secure?
    start = int(request.args['start'])
    length = int(request.args['length'])

    # cs is custom search
    csKey = request.args['csKey']
    csField = request.args['csOpt']
    retlist = []
    res_csv = pd.read_csv("%s%s.csv"%(app.config['UPLOAD_FOLDER'],task_id)).values.tolist()
    for i in range(start,start+length):
        row = res_csv[i]
        if not filter_search(row,csKey,csField):
            continue
        wild = row[1][:5] + '<span class="bolded-red">' + row[1][5] + '</span>' + row[1][6:]
        mut = row[2][:5] + '<span class="bolded-red">' + row[2][5] + '</span>' + row[2][6:]
        p_or_z_score = "%.3e"%row[4] if abs(row[4]) < 0.0001 else "%.4f"%row[4]
        retlist.append([int(row[0]),wild,mut,"%.4f"%row[3],p_or_z_score] + row[5:])

    # check orderable -- we disable orderMulti in result.js so we can assume
    # one column ordering.
    order_col = int(request.args["order[0][column]"])
    order_reverse = False if request.args["order[0][dir]"] == "asc" else True
    retlist = sorted(retlist, key = lambda x: x[order_col],reverse=order_reverse)

    return jsonify({
        "draw": draw,
        "recordsTotal": len(retlist),
        "recordsFiltered": len(retlist),
        "data": retlist
    })


'''
task_id: specifies unique url for each chain, this is the same with the id of
         the youngest child in the pipeline
in javascript: updateProgress function calls this using status url
'''
@app.route('/status/<task_id>')
def task_status(task_id):
    pretask_id = request.args["parent-0"] # i.e. preprocess

    # see which task is active:
    task = celery.AsyncResult(pretask_id) # preprocess
    if task.state == 'SUCCESS': # if preprocess is finished then check prediction
        task = celery.AsyncResult(task_id) # prediction
    if task.state == 'PENDING':
        response = {
            'state': task.state,
            'current': 0,
            'total': 1,
            'status': 'Pending...'
        }
        task.forget()
    else: #task.state == 'PROGRESS' or 'SUCCESS'?
        response = {
            'state': task.state,
            'current': task.info.get('current', 0),
            'total': task.info.get('total', 1),
            'status': task.info.get('status', '')
        }
        if 'result' in task.info:
            response['result'] = task.info['result']
            response['csvlink'] = task.info['csvlink']
        # TODO: need to forget task

    #task.forget() #??? not sure if needed
    # pkill -9 -f 'celery worker'

    return jsonify(response)
    #return render_template("result.html",tblres=[]) #tblres=res
    #return send_from_directory(app.config['UPLOAD_FOLDER'],
    #                           filename)

@app.route('/getinputparam/<job_id>', methods=['GET'])
def get_input_param(job_id):
    # key: filename,pbmselected,filteropt,filterval,chrver
    indict = db.hgetall(job_id)
    indict["genes_selected"] = ast.literal_eval(indict["genes_selected"])
    return json.dumps(indict) # must return a json

# ==========================

'''
return empty string if there is no file, filename if there is file and it has been handled
by the function
return: status, msg
msg = filename if success
'''
def prepare_request(request):
    # There is no input file in the request
    if 'input-file' not in request.files:
        return 'error','no input file part'
    file = request.files['input-file']
    # Can only accept tsv or csv
    if not allowed_file(file.filename):
        return 'error','please upload only tsv or csv'
    # No file selected:
    if file.filename == '':
        return 'error','no selected file'
    # No TFs selected
    if not request.form.getlist('pred-select'):
        return 'error','please select transcription factors'
    # Check if we have all the columns we need
    filename = secure_filename(file.filename)
    filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
    file.save(filepath)
    required = ['chromosome','chromosome_start','mutation_type','mutated_from_allele','mutated_to_allele']
    with open(filepath) as f:
        file_extension = os.path.splitext(filepath)[1]
        if file_extension == ".tsv":
            incols = f.readline().strip().split("\t")
        else: # must be csv since we checked it
            incols = f.readline().strip().split(",")
    if not all(elem in incols for elem in required):
        os.remove(filepath)
        return 'error','some required fields are missing from the input file'
    # Finally, file is okay
    return 'success',filename


@app.route('/upload', methods=['POST'])
def handle_upload():
    tfpref = "prediction6mer." # for rapid, need to be empty for now
    tfext = ".csv"
    if request.method == 'POST':
        status,msg = prepare_request(request)
        if status=='error':
            return jsonify({'Message':msg}), 500
        else:
            # request.form.getlist('pred-select'):['Arid3a:Arid3a_3875.1_v1_deBruijn', 'Bhlhb2:Bhlhb2_4971.1_v1_deBruijn']
            genes_selected = [elm.split(":")[0] for elm in request.form.getlist('pred-select')]

            select_list = [elm.split(":")[1] for elm in request.form.getlist('pred-select')]
            unique_pbms = list({tfpref+x+tfext for pbm in select_list for x in pbm.split(',')})

            chrver = request.form.get('genome-select')

            filteropt = int(request.form.get('optradio'))
            if filteropt == 1:
                filterval = int(request.form.get('output-selection-opt'))
            else:
                filterval = float(request.form.get('output-selection-opt'))

            task = chain(celerytask.inittbl.s(app.config['UPLOAD_FOLDER'] + msg,
                        app.config['CHRDIR'] +"/"+chrver),
                        celerytask.do_prediction.s(unique_pbms,genes_selected,filteropt,filterval)).apply_async() # put genes_selected here

            # better to not use nested dict as Flask use CallbackDict to track modification
            session['%s_p0'%task.id] = task.parent.id
            session['%s_p1'%task.id] = task.id

            # ==== STORING IN REDIS PART ====
            session_info = {"filename":msg,
                            "genes_selected":genes_selected,
                            "filteropt":filteropt,
                            "filterval":filterval,
                            "chrver":request.form.get('genome-select')}
            if db.exists(task.id):
                db.delete(task.id)
            db.hmset(task.id,session_info)
            db.expire(task.id, app.config['USER_DATA_EXPIRY'])
            # ================================

            task.forget() # not sure if needed???

            resp = make_response(jsonify({}), 202, {'Location': url_for('process_request',job_id=task.id)})

            job_name = request.form.get("job-name") if request.form.get("job-name") else task.id
            resp.set_cookie("qbic_recents:%s"%task.id, job_name, max_age=app.config['USER_DATA_EXPIRY'])
            return  resp # {'Location': url_for('task_status',task_id=task.id)

            #return redirect(url_for('process_request'),code=202)'''

#================================

@app.route('/predlist', methods=['GET'])
def get_predlist():
    with open(app.config['HUGO_PBM_MAPPING'],'r') as f:
        family_map = {}
        for line in f:
            key,val = line.strip().split("->")
            valmap = {z[0]:z[1] for z in (y.split(":") for y in (x for x in val.split(";")))} # generator
            family_map[key] = valmap
    return jsonify(family_map)

@app.route('/recent', methods=['GET'])
def get_recent_jobs():
    recents = [key for key in request.cookies.keys() if key.startswith("qbic_recents:")]
    if recents:
        rj_urls = [[request.cookies.get(job_key),url_for('process_request',job_id=job_key.split(":",1)[1])] for job_key in recents]
        return json.dumps(rj_urls)
    else:
        return json.dumps({})

#==============================================


@app.route('/', methods=['GET', 'POST'])
def index():
    #session.permanent = True
    #session.clear() -- need to limit the amount of session somewhere
    return render_template("index.html")

'''
job_id: from last task in the pipeline
'''
@app.route('/process/<job_id>', methods=['GET'])
def process_request(job_id):
    # get the information saved by handle_upload
    p0 = session['%s_p0'%job_id]
    p1 = session['%s_p1'%job_id]
    parents = json.dumps({'parent-0':p0,'parent-1':p1})
    #session.clear() # clear the session given from index
    return render_template("result.html",stats_url=url_for('task_status',task_id=job_id),parents=parents)

@app.route('/about')
def about():
    return render_template("about.html")
