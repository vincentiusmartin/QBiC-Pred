from flask import render_template
from app import app

@app.route('/', methods=['GET', 'POST'])
def index():
    #session.permanent = True
    #session.clear() -- need to limit the amount of session somewhere
    return render_template("index.html")

@app.route('/about')
def about():
    return render_template("about.html")

@app.route('/downloads')
def downloads():
    return render_template("downloads.html")
