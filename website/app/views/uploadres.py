

@app.route('/makepred', methods=['GET', 'POST'])
def makepred():
    #session.permanent = True
    #session.clear() -- need to limit the amount of session somewhere
    return render_template("makepred.html")
