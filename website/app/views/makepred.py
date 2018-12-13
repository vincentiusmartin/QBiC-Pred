import os
import sys
sys.path.insert(0, '..')

from flask import request,render_template,make_response,jsonify,url_for
from werkzeug.utils import secure_filename

from celery import Celery,chain
from app import app,db

import app.controller.celerytask as celerytask

ALLOWED_EXTENSIONS = set(['csv','tsv'])
#app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

# view specific utils
def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS
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
    # Check file size
    file.seek(0, os.SEEK_END)
    file_length = file.tell()
    if file_length > app.config['MAX_FILE_LENGTH']:
        maxsize = app.config['MAX_FILE_LENGTH'] / (1024*1024)
        return 'error','uploaded file is larger than the allowed maximum size of %dMB' % maxsize
    file.seek(0) # seek back
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
    tfext = ".txt"
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

            # ==== STORING IN REDIS PART ====
            # it is important to store these in redis so information can be
            # passed to different browsers/machines.
            session_info = {"parent_id":task.parent.id,
                            "task_id":task.id,
                            "filename":msg,
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
            # we can put this in cookie to let browser save the recent jobs
            resp.set_cookie("qbic_recents:%s"%task.id, job_name, max_age=app.config['USER_DATA_EXPIRY'])
            return resp # {'Location': url_for('task_status',task_id=task.id)

            #return redirect(url_for('process_request'),code=202)'''

@app.route('/makepred', methods=['GET', 'POST'])
def makepred():
    #session.permanent = True
    #session.clear() -- need to limit the amount of session somewhere
    return render_template("makepred.html")

@app.route('/about')
def about():
    return render_template("about.html")
