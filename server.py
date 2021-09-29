import os

from flask import Flask, render_template, flash, request, redirect, url_for
from werkzeug.utils import secure_filename

UPLOAD_FOLDER = 'data/'
ALLOWED_EXTENSIONS = {'xlsx'}

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
import step1

@app.route('/')
def index():
  return render_template('index.html')

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@app.route('/', methods=['POST'])
def upload_file():
  if request.method == 'POST':
    # check if the post request has the file part
    if 'file' not in request.files:
      flash('No file part')
      return redirect(request.url)
    file = request.files['file']
    # If the user does not select a file, the browser submits an
    # empty file without a filename.
    if file.filename == '':
      flash('No selected file')
      return redirect(request.url)
    if file and allowed_file(file.filename):
      filename = secure_filename(file.filename)
      file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
      platedf = step1.delimit_insertion(os.path.join(app.config['UPLOAD_FOLDER'], filename))
      platedf.to_csv(os.path.join(app.config['UPLOAD_FOLDER'],filename.split(".xlsx")[0]+'primers_in.csv'), columns=['sample', 'bed_range'], header=None, index=None)
      return file.filename
      # return redirect(url_for('_file', name=filename))
  # return

# @app.route('/my-link/')
# def my_link():
#   step1.delimit_insertion()

if __name__ == '__main__':
  app.run(debug=True)