import sys
sys.path.insert(0, '..')

from flask import send_from_directory,jsonify,request,render_template,url_for,Response,abort
from app import app,db,celery

import redisearch

# ast is used to convert string literal representation of list to a list,
# this is needed since redis stores list as string
import json,ast

import pandas as pd

'''
job_id: from last task in the pipeline
'''
@app.route('/process/<job_id>', methods=['GET'])
def process_request(job_id):
    # get the task data from redis
    taskdata = db.hgetall(job_id)
    if "parent_id" not in taskdata:
        abort(404)
    p0 = taskdata["parent_id"]
    p1 = taskdata["task_id"]
    parents = json.dumps({'parent-0':p0,'parent-1':p1})
    #session.clear() # clear the session given from index
    return render_template("result.html",stats_url=url_for('task_status',task_id=job_id),parents=parents)

@app.route('/files/<taskid>')
def get_file(taskid):
    """Download a file."""
    task = db.hgetall("%s:cols"%taskid)
    cols = ast.literal_eval(task['cols'])
    csv = "%s\n" % ",".join(cols)

    client = redisearch.Client(taskid)
    num_docs = int(client.info()['num_docs'])

    for i in range(0,num_docs):
        doc = client.load_document(i)
        listparam = []
        for col in cols:
            if col == "TF_gene":
                listparam.append("\"%s\""eval("doc.%s"%col))
            else:
                listparam.append(eval("doc.%s"%col))
        csv += ",".join(listparam)
        if i != num_docs-1:
            csv += "\n"

    ''' return the csv file without having to save it '''
    return Response(
        csv,
        mimetype="text/csv",
        headers={"Content-disposition":
                 "attachment; filename=prediction_result-%s.csv"%taskid})

@app.route('/getrescol/<task_id>',methods=['GET'])
def get_res_col(task_id):
    cols = []
    col_id = "%s:cols"%task_id
    if db.exists(col_id):
        cols_fromdb = ast.literal_eval(db.hgetall(col_id)['cols'])
        #with open("%s%s.csv"%(app.config['UPLOAD_FOLDER'],task_id)) as f:
        cols = [{"title":title} for title in cols_fromdb]
    # else just return an empty list
    return jsonify(cols)

def dofilter(search_filter,doc):
    flag = True
    for filter in search_filter:
        if filter["searchOpt"] == "in sequence":
            if filter["searchKey"] in doc.wild or filter["searchKey"] in doc.mutant:
                continue
        elif filter["searchOpt"] == "at least" or filter["searchOpt"] == "at most":
            col = filter["searchCol"].replace("-","_")
            searchval = float(getattr(doc,col))
            threshold = float(filter["searchKey"])
            if (filter["searchOpt"] == "at least" and searchval >= threshold) or \
                (filter["searchOpt"] == "at most" and searchval <= threshold):
                continue
        else:
            searchval = getattr(doc,filter["searchCol"])
            if (filter["searchOpt"] == "exact" and filter["searchKey"] == searchval) or \
               (filter["searchOpt"] == "exclude" and filter["searchKey"] != searchval):
                continue
        flag = False
        break
    return flag

def filter_fromdb(task_id,search_filter,start,length,order_col="row",order_asc=True):
    result = {}
    client = redisearch.Client(task_id)
    result['recordsTotal'] = int(client.info()['num_docs'])

    #manual = False # manually made because redisearch sucks
    if search_filter:
        query = redisearch.Query("*").sort_by(order_col,order_asc).paging(0,result['recordsTotal'])
        documents = client.search(query).docs
        filtered_docs = list(filter(lambda doc: dofilter(search_filter,doc),documents))
        result['recordsFiltered'] = len(filtered_docs)
        result['data'] = filtered_docs[start:start+length]
    else:
        query = redisearch.Query("*").sort_by(order_col,order_asc).paging(start,length)
        res = client.search(query)
        result['recordsFiltered'] = res.total
        result['data'] = res.docs

    #if searchtype == "exclude":
    #    #searchtext = searchquery.replace(">","\\>") -- @col:query
    #    querystr = "-(@%s:\"%s\")" % (colname,searchquery)
    #else: #searchtype == "exact"
    #    querystr = "@%s:\"%s\"" % (colname,searchquery)

    return result

def customround(num):
    fl = float(num) # so it can accept string
    return "%.3e"%fl if abs(fl) < 10**(-4) or abs(fl) > 10**(4) else "%.4f"%fl

def htmlformat(invar,type,colname):
    str_in = str(invar)
    if type == "filter":
        return """\
          <div class="dropdown cell-filter">
            <button class="btn btn-link unstyled-button cell-btn" type="button" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
              {content}
            </button>
            <div class="dropdown-menu" aria-labelledby="dropdownMenuButton">
              <button class="dropdown-item cell-filter-item" data-colname={colname} data-filter="exact">Only include rows with this value</a>
              <button class="dropdown-item cell-filter-item" data-colname={colname} data-filter="exclude">Exclude rows with this value</a>
            </div>
          </div>
        """.format(content=str_in,colname=colname)

@app.route('/getrestbl/<task_id>',methods=['GET'])
def get_res_tbl(task_id):
    # Get all mandatory informations we need
    # Info about params: https://datatables.net/manual/server-side#DataTables_Table_1
    draw = int(request.args['draw']) # not secure # TODO: make it secure?
    start = int(request.args['start'])
    length = int(request.args['length'])

    # cs is custom search
    searchFilter = ast.literal_eval(request.args['searchFilter'])

    # check orderable -- orderMulti is disabled in result.js so we can assume
    # one column ordering.
    order_col = int(request.args["order[0][column]"])
    order_asc = True if request.args["order[0][dir]"] == "asc" else False

    retlist = []
    #res_csv = pd.read_csv("%s%s.csv"%(app.config['UPLOAD_FOLDER'],task_id)).values.tolist()

    cols = ast.literal_eval(db.hgetall("%s:cols"%task_id)['cols'])

    filtered_db = filter_fromdb(task_id,searchFilter,start,length,cols[order_col],order_asc)
    retlist = []

    for doc in filtered_db['data']:
        #if not filter_search(row,csKey,csField):
        #    continue
        rowdict = {col:getattr(doc,col) for col in cols}
        rowdict['wild'] = doc.wild[:5] + '<span class="bolded-red">' + doc.wild[5] + '</span>' + doc.wild[6:]
        rowdict['mutant'] = doc.mutant[:5] + '<span class="bolded-red">' + doc.mutant[5] + '</span>' + doc.mutant[6:]
        rowdict['z_score'] = customround(doc.z_score)
        rowdict['p_value'] = customround(doc.p_value)
        rowdict['binding_status'] = htmlformat(doc.binding_status,"filter","binding_status")
        rowdict['gapmodel'] = htmlformat(doc.gapmodel,"filter","gapmodel")
        rowdict['TF_gene'] = htmlformat(doc.TF_gene,"filter","TF_gene")
        rowdict['pbmname'] = htmlformat(doc.pbmname,"filter","pbmname")
        retlist.append([rowdict[col] for col in cols])

    return jsonify({
        "draw": draw,
        "recordsTotal": filtered_db['recordsTotal'],
        "recordsFiltered": filtered_db['recordsFiltered'],
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
    response = {}
    if pretask_id == "uploadpred":
        response = {
            'state': 'SUCCESS',
            'current': 1,
            'total': 1,
            'status': '',
            'result':'',
            'taskid':task_id
        }
    else:
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
                response['taskid'] = task.info['taskid']
            # TODO: need to forget task

        #task.forget() #??? not sure if needed
        # pkill -9 -f 'celery worker'
    return jsonify(response)

@app.route('/getinputparam/<job_id>', methods=['GET'])
def get_input_param(job_id):
    # key: filename,pbmselected,filteropt,filterval,chrver
    indict = db.hgetall(job_id)
    indict["genes_selected"] = ast.literal_eval(indict["genes_selected"])
    return json.dumps(indict) # must return a json
