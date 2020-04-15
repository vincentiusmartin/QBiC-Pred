import sys
sys.path.insert(0, '..')

from flask import send_from_directory,jsonify,request,render_template,url_for,Response,abort
from app import app,redisdb,celery,mongodb

import redisearch
import billiard # just for contrive error catching

# ast is used to convert string literal representation of list to a list,
# this is needed since redis stores list as string
import json,ast
import pandas as pd
import re

'''
job_id: from last task in the pipeline
'''
@app.route('/process/<job_id>', methods=['GET'])
def process_request(job_id):
    # get the task data from redis
    taskdata = redisdb.hgetall(job_id)
    if "parent_id" not in taskdata:
        abort(404)
    p0 = taskdata["parent_id"]
    p1 = taskdata["task_id"]
    parents = json.dumps({'parent-0':p0,'parent-1':p1})
    #session.clear() # clear the session given from index
    return render_template("result.html",stats_url=url_for('task_status',task_id=job_id),parents=parents)

# /<taskid>/<filters>
@app.route('/files/<filetype>/<task_id>/<filters>')
def get_file_fromtbl(filetype,task_id,filters): #taskid,filters
    #filtered_db = filter_fromdb(task_id,searchFilter,start,length,cols[order_col],order_asc)
    search_filter = ast.literal_eval(filters); # [{"searchOpt":"in sequence","searchKey":"AAT","searchCol":"sequence"}
    filtered = filter_fromdb(task_id,search_filter,start=0,length=-1)

    ftype = filetype.lower()
    if ftype == "tsv":
        sep = "\t"
    else:
        sep = ","

    cols = get_mongocols(task_id)
    tblret = sep.join(cols) + "\n"
    for doc in filtered['data']:
        try: # TODO: handle this
            row = []
            for col in cols: # need to do this because we can have multiple TF_genes
                if col == "TF_gene":
                    row.append("\"%s\""%doc[col])
                elif col == "p_value":
                    row.append('{:0.3e}'.format(doc[col]))
                else:
                    row.append(str(doc[col]))
            tblret += sep.join(row) + "\n"
        except Exception as e: #now: if not found, just return 404
            print("Exception: " + str(e))
            abort(404)
    return Response(
        tblret[:-1],
        mimetype="text/csv",
        headers={"Content-disposition":
                 "attachment; filename=prediction_result-%s.%s"%(task_id,ftype)})

@app.route('/filesdb/<filetype>/<taskid>')
def get_file_fromdb(filetype,taskid):
    """Download a file."""
    #task = db.hgetall("%s:cols"%taskid)
    cols = get_mongocols(taskid)
    sep = "\t" if filetype=="tsv" else ","
    csv = "%s\n" % sep.join(cols)

    client = redisearch.Client(taskid)
    num_docs = int(client.info()['num_docs'])

    for i in range(0,num_docs):
        doc = client.load_document("%s_%d"%(taskid,i))
        listparam = []
        for col in cols:
            if col == "TF_gene":
                listparam.append("\"%s\""%eval("doc.%s"%col))
            else:
                listparam.append(eval("doc.%s"%col))
        csv += sep.join(listparam)
        if i != num_docs-1:
            csv += "\n"

    ''' return the csv/tsv file without having to save it '''
    return Response(
        csv,
        mimetype="text/csv",
        headers={"Content-disposition":
                 "attachment; filename=prediction_result-%s.%s"%(taskid,filetype)})

def get_mongocols(task_id):
    if task_id in mongodb.list_collection_names():
        # collection = mongodb[task_id]
        # cols = list(collection.find_one().keys()) # get column name
        # cols.remove("_id")
        # return cols
        return ['row', 'wild', 'mutant', 'diff', 'z_score', 'p_value', 'binding_status', 'TF_gene', 'pbmname']
    else: # just return an empty list
        return []

@app.route('/getrescol/<task_id>',methods=['GET'])
def get_res_col(task_id):
    colmap = { # hardcoded, need to fix this
    "row":"Index",
    "wild":"Ref",
    "mutant":"Alt",
    "diff":"Difference",
    "p_value":"p value",
    "z_score":"z score",
    "TF_gene":"TF gene",
    "binding_status":"Binding status",
    #"gapmodel":"Gap model", #vmartin: comment this
    "pbmname":"PBM filename"
    }
    cols_fromdb =  get_mongocols(task_id)
    cols = []
    orderable_cols = ['diff','row','z_score','p_value']
    for title in cols_fromdb:
        # we only provide sort on some columns
         #vmartin: hide gap model
        if title in colmap:
            if title not in orderable_cols:
                cols.append({"title":colmap[title], "orderable": 0})
            else:
                cols.append({"title":colmap[title], "orderable": 1})
    return jsonify(cols)

def query_filter(search_filter):
    query_or = {}
    query = {}
    inseq_substr = ""
    exact_cols = []
    for q in search_filter:
        print(q["searchOpt"])
        if q["searchOpt"] == "in sequence":
            inseq_substr += "%s|" % q["searchKey"]
        elif  q["searchOpt"] == "or":
            if q["searchCol"] not in query_or:
                query_or[q["searchCol"]] = []
            query_or[q["searchCol"]].append({q["searchCol"]:q["searchKey"]})
        elif q["searchOpt"] == "at most" or q["searchOpt"] == "at least":
            op = "$lte" if  q["searchOpt"] == "at most" else "$gte"
            thres = float(q["searchKey"])
            if q["searchCol"] == "p-value":
                query["p_value"] = {op:thres}
            elif q["searchCol"] == "z-score":
                # use absolute value if z-score, for now we only do this for z_score, so
                # I think it's fine to use expr here
                query["$expr"] =  {op: [ {"$abs": "$z_score"} , abs(thres) ] }
        elif q["searchOpt"] == "exact":
            if q["searchCol"] not in query:
                exact_cols.append(q["searchCol"])
                query[q["searchCol"]] = ""
            query[q["searchCol"]] += "%s|" % q["searchKey"]

    # change all exact to regex
    for col in exact_cols:
        query[col] = {"$regex":re.compile(query[col][:-1], re.I)}

    # make the query
    query_and = [{"$or":v} for k,v in query_or.items() if k != "z-score" and k != "p-value"]

    # take care of the string sequence
    if inseq_substr: # if the query exist
        inseq_substr = inseq_substr[:-1]
        pat = re.compile(inseq_substr, re.I)
        query_str = {"$or":[{"wild":{"$regex":inseq_substr}}, {"mutant":{"$regex":inseq_substr}}]}
        query_and.append(query_str)

    if query_and:
        query["$and"] = query_and
    return query

def filter_fromdb(task_id,search_filter,start,length=-1,order_col="row",order_asc=1):
    '''
    task_id,
    search_filter,
    start,length = -1,
    order_col="row",
    order_asc=True
    '''
    result = {}
    collection = mongodb[task_id]
    result['recordsTotal'] = collection.count()

    #manual = False # manually made because redisearch sucks
    if length == -1:
        length = result['recordsTotal'] - start

    # if there is filter or length == -1 we return everything
    # hay que devolver todo porque necesitamos contar el número de filas
    if search_filter:
        #filtered_docs = collection.find({}, {'_id': False}).sort(order_col,order_asc)
        #filtered = dofilter(search_filter,collection)
        query = query_filter(search_filter) # count_documents()
        documents = collection.find(query, {'_id': False}) #.sort(order_col,order_asc)
        result['recordsFiltered'] = documents.count()
        sorted_res = documents.sort(order_col,order_asc) \
                               .skip(start) \
                               .limit(length)
    else: # list(collection.find({}, {'_id': False}))
        # if there is no filter, get the paging
        result['recordsFiltered'] = result['recordsTotal']
        sorted_res = collection.find({}, {'_id': False}) \
                               .sort(order_col,order_asc) \
                               .skip(start) \
                               .limit(length)
    result['data'] = list(sorted_res)

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
    '''
    invar is the value of the cell
    '''
    str_in = str(invar)
    if type == "filter":
        buttonhtml = """\
            <span class="dropdown cell-filter">
              <button class="btn btn-link unstyled-button cell-btn" type="button" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
                {content}
              </button>
              <div class="dropdown-menu" aria-labelledby="dropdownMenuButton">
                <button class="dropdown-item cell-filter-item" data-colname={colname} data-filter="exact">Only include rows with this value</button>
                <button class="dropdown-item cell-filter-item" data-colname={colname} data-filter="exclude">Exclude rows with this value</button>
                {additional}
              </div>
            </span>
        """
        buttonhgnc = """<a href="https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{hgnc}" target="_blank"  class="btn btn-link dropdown-item" role="button">See gene entry on HGNC website</a>"""
        buttonmodels = """<a href="/models?search={search}" target="_blank"  class="btn btn-link dropdown-item" role="button">Show OLS models for this gene</a>"""
        buttonsearchpbm = """<a href="/models?search={search}" target="_blank"  class="btn btn-link dropdown-item" role="button">Show model information</a>"""
        buttondownloadpbm =  """<a href="/download/pbmdata/{pbmname}"  class="btn btn-link dropdown-item" role="button">Download PBM data</a>"""

        if colname == "TF_gene":
            cellvals = str_in.split(",")
            formatted = []
            for val in cellvals:
                # if we don't have the entry for the gene, then just fill empty
                # This is useful when we have gene that is used for testing only
                if val in app.config['HUGO_NAME_ID_MAPPING']:
                    hgncstrbtn = buttonhgnc.format(hgnc=app.config['HUGO_NAME_ID_MAPPING'][val])
                    modelstrbtn = buttonmodels.format(search=val)
                else:
                    hgncstrbtn = buttonhgnc.format(hgnc="")
                    modelstrbtn = buttonmodels.format(search="")
                btnetc = hgncstrbtn + modelstrbtn
                formatted.append(buttonhtml.format(content=val,colname=colname,additional=btnetc))
            content = ""
            charinrow = 0
            total = len(formatted)
            for i in range(0,total): #super forced hacking for string wrap -_-
                content += formatted[i]
                if i < total-1:
                    content += ", "
                charinrow += len(cellvals[i]) + 2
                if charinrow > 20 and i<total-1:
                    content += "<br />"
                    charinrow = 0
        elif colname == "pbmname":
            pbmstrbutton = buttondownloadpbm.format(pbmname="%s.txt"%str_in)
            modelstrbtn = buttonsearchpbm.format(search=str_in)
            btnetc = pbmstrbutton + modelstrbtn
            content = buttonhtml.format(content=str_in,colname=colname,additional=btnetc)
        else:
            content = buttonhtml.format(content=str_in,colname=colname,additional="")

        return content

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
    if request.args['columns[%d][orderable]' % order_col] == '0':
        order_col = 0
        order_asc = 1
    else:
        order_asc = 1 if request.args["order[0][dir]"] == "asc" else -1

    retlist = []
    #res_csv = pd.read_csv("%s%s.csv"%(app.config['UPLOAD_FOLDER'],task_id)).values.tolist()


    cols = get_mongocols(task_id) # ast.literal_eval(

    filtered_db = filter_fromdb(task_id,searchFilter,start,length,cols[order_col],order_asc)
    # filtered = filter_fromdb(task_id,search_filter,start=0,length=-1) -- in download

    retlist = []

    # for now, if nothing found then return an empty table
    if filtered_db['data'] and len(filtered_db['data']) == 0: #not hasattr(filtered_db['data'][0],'row'): # hasattr
        return jsonify({
            "draw": 0,
            "recordsTotal": 0,
            "recordsFiltered": 0,
            "data": [],
            "status":"false"
        })

    for doc in filtered_db['data']:
        #if not filter_search(row,csKey,csField):
        #    continue
        #try: # TODO: handle this
        #    rowdict = {col:getattr(doc,col) for col in cols}
        #except: #now: if not found, just return 404
        #    abort(404)
        rowdict = {}
        rowdict['row'] = doc['row']
        rowdict['wild'] = doc['wild'][:5] + '<span class="bolded-red">' + doc['wild'][5] + '</span>' + doc['wild'][6:]
        rowdict['mutant'] = doc['mutant'][:5] + '<span class="bolded-red">' + doc['mutant'][5] + '</span>' + doc['mutant'][6:]
        rowdict['diff'] = doc['diff']
        rowdict['z_score'] = customround(doc['z_score'])
        rowdict['p_value'] = customround(doc['p_value'])
        rowdict['binding_status'] = htmlformat(doc['binding_status'],"filter","binding_status") #vmartin: binding-flag
        #rowdict['gapmodel'] = htmlformat(doc.gapmodel,"filter","gapmodel")
        rowdict['TF_gene'] = htmlformat(doc['TF_gene'],"filter","TF_gene")
        rowdict['pbmname'] = htmlformat(doc['pbmname'],"filter","pbmname")
        retlist.append([rowdict[col] for col in cols])

    print(cols)
    return jsonify({
        "draw": draw,
        "recordsTotal": filtered_db['recordsTotal'],
        "recordsFiltered": filtered_db['recordsFiltered'],
        "data": retlist,
        "status":"success"
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
            # to handle worker die
            # TODO: detect when worker dies and handle
            if isinstance(task.info, billiard.exceptions.WorkerLostError):
                print("Worker is killed for %s, returning an error message to the user..." % str(task_id))
                response = {'state': 'ERROR',
                            'error': 'There was an error running the job. The webserver authors have recorded the error and will work to address it. Please re-submit your job. We apologize for the inconvenience!'
                            }
            else:
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
    return jsonify(response)

@app.route('/getinputparam/<job_id>', methods=['GET'])
def get_input_param(job_id):
    # key: filename,pbmselected,filteropt,filterval,chrver
    indict = redisdb.hgetall(job_id)
    indict["genes_selected"] = ast.literal_eval(indict["genes_selected"])
    return json.dumps(indict) # must return a json
