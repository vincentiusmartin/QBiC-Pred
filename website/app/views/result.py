import sys
sys.path.insert(0, '..')

from flask import send_from_directory,jsonify,request,render_template,url_for
from app import app,db,celery

import json,ast

import pandas as pd

def filter_search(infields,keyword,type):
    if type == "All":
        if any(keyword in str(field) for field in infields):
            return True
        else:
            return False
    return True

'''
job_id: from last task in the pipeline
'''
@app.route('/process/<job_id>', methods=['GET'])
def process_request(job_id):
    # get the task data from redis
    taskdata = db.hgetall(job_id)
    p0 = taskdata["parent_id"]
    p1 = taskdata["task_id"]
    parents = json.dumps({'parent-0':p0,'parent-1':p1})
    #session.clear() # clear the session given from index
    return render_template("result.html",stats_url=url_for('task_status',task_id=job_id),parents=parents)

@app.route('/files/<path:filename>')
def get_file(filename):
    """Download a file."""
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename,as_attachment=True)

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

    filtered_tbl = [row for row in res_csv if filter_search(row,csKey,csField)]

    for i in range(start,min(len(filtered_tbl),start+length)):
        row = filtered_tbl[i]
        #if not filter_search(row,csKey,csField):
        #    continue
        wild = row[1][:5] + '<span class="bolded-red">' + row[1][5] + '</span>' + row[1][6:]
        mut = row[2][:5] + '<span class="bolded-red">' + row[2][5] + '</span>' + row[2][6:]
        p_or_z_score = "%.3e"%row[4] if abs(row[4]) < 10**(-4) and abs(row[4]) > 10**(-10) else "%.4f"%row[4]
        # str(elm) for elm in row[5:] -> used to fix nan
        retlist.append([int(row[0]),wild,mut,"%.4f"%row[3],p_or_z_score] + [str(elm) for elm in row[5:]])
    # check orderable -- we disable orderMulti in result.js so we can assume
    # one column ordering.
    order_col = int(request.args["order[0][column]"])
    order_reverse = False if request.args["order[0][dir]"] == "asc" else True
    retlist = sorted(retlist, key = lambda x: x[order_col],reverse=order_reverse)

    print(retlist)

    return jsonify({
        "draw": draw,
        "recordsTotal": len(res_csv),
        "recordsFiltered": len(filtered_tbl),
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
        if 'error' in task.info:
            response['error'] = task.info['error']
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
