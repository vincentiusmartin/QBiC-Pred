from flask import render_template
from app import app

@app.route('/', methods=['GET', 'POST'])
def index():
    #session.permanent = True
    #session.clear() -- need to limit the amount of session somewhere
    return render_template("index.html")
