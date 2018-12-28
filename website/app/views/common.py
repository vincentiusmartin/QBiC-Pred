from flask import send_from_directory

from app import app

@app.route('/download/<path:filename>')
def download_example(filename):
    return send_from_directory(app.config['STATIC_EXAMPLE_DIR'], filename, as_attachment=True)
