from flask import render_template,jsonify
from app import app

import pandas as pd

@app.route('/models')
def models():
    return render_template("models.html")

#@app.route('/getmodelscol',methods=['GET'])
#def get_models_col():
#    cols = pd.read_csv(app.config['MODELS_TBL_PATH'], nrows=0).columns.tolist()
    #df = pd.read_csv(app.config['MODELS_TBL_PATH'],dtype=str).drop(['best'], axis=1)
    #df.to_csv("test.csv",index=False)
    # we need "title to indicate it is a column"
#    return jsonify([{"title":col} for col in cols])

@app.route('/getmodeltbl',methods=['GET','POST'])
def get_models_tbl():
    # data frame cannot have na...
    df = pd.read_csv(app.config['MODELS_TBL_PATH'],dtype=str,keep_default_na=False)
    df = df.rename(columns={col: col.replace("."," ").replace("_"," ") for col in df.columns})
    #df = df.head(n=200)
    dict_df = df.to_dict('records')
    # return should be in data object to be parsed by DataTables
    jsonified = jsonify({
        "cols":list(df.columns.values), # to preserve the order of the colnames
        "data":dict_df
    })
    return jsonified
