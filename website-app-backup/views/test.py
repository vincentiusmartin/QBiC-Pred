
from flask import *

from app import app,celery

# for testing only
@app.route('/testing', methods=['GET'])
def testing():
    return render_template("test.html") # ,test_url=url_for('get_test_tbl')

@app.route('/gettestcol',methods=['GET'])
def test_get_col():
    with open("/Users/vincentiusmartin/Research/MutationPredictor/tfbc-website/ed19bc67-e5b1-46bd-b250-9efd3a7cc972.csv") as f:
        cols = f.readline().strip().split(",")
    return jsonify(cols)

@app.route('/gettesttbl',methods=['GET'])
def test_get_tbl():
    # {\"0\":\"rowidx\",\"1\":\"wild-type\",\"2\":\"mutant\",\"3\":\"diff\",\"4\":\"z-score\",\"5\":\"pbmname\"},
    draw = int(request.args['draw']) # not secure # TODO: make it secure
    start = int(request.args['start'])
    length = int(request.args['length'])
    retlist = []
    with open("/Users/vincentiusmartin/Research/MutationPredictor/tfbc-website/ed19bc67-e5b1-46bd-b250-9efd3a7cc972.csv") as f:
        next(f)
        for line in f:
            splitted = line.strip().split(",")
            row = [splitted[i] for i in range(0,len(splitted))]
            retlist.append(row)
    return jsonify({
        "draw": draw,
        "recordsTotal": len(retlist),
        "recordsFiltered": len(retlist),
        "data": retlist[start:start+length]
    })
    '''return jsonify({
  "draw": 1,
  "recordsTotal": 57,
  "recordsFiltered": 57,
  "data": [
    [
      "Airi",
      "Satou",
      "Accountant",
      "Tokyo",
      "28th Nov 08",
      "$162,700"
    ],
    [
      "Angelica",
      "Ramos",
      "Chief Executive Officer (CEO)",
      "London",
      "9th Oct 09",
      "$1,200,000"
    ],
    [
      "Ashton",
      "Cox",
      "Junior Technical Author",
      "San Francisco",
      "12th Jan 09",
      "$86,000"
    ],
    [
      "Bradley",
      "Greer",
      "Software Engineer",
      "London",
      "13th Oct 12",
      "$132,000"
    ],
    [
      "Brenden",
      "Wagner",
      "Software Engineer",
      "San Francisco",
      "7th Jun 11",
      "$206,850"
    ],
    [
      "Brielle",
      "Williamson",
      "Integration Specialist",
      "New York",
      "2nd Dec 12",
      "$372,000"
    ]]})'''
