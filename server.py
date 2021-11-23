from datetime import datetime
import os
import tempfile
from os.path import exists

import boto3
import paramiko
from flask import Flask
from flask import render_template, flash, request, redirect
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
        if 'xlsxfile' not in request.files:
            print('No file part')
            return redirect(request.url)
        file = request.files['xlsxfile']
        # If the user does not select a file, the browser submits an
        # empty file without a filename.
        if file.filename == '':
            print('No selected file')
            return redirect(request.url)
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            #     file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            #     platedf = step1.delimit_insertion(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            app.config['CSV'] = os.path.join(app.config['UPLOAD_FOLDER'], filename.split(".xlsx")[0] + 'primers_in.csv')
        #     platedf.to_csv(app.config['CSV'], columns=['sample', 'bed_range'], header=None, index=None)

        if exists(app.config['CSV']):
            flash("Coordinates csv file successfully generated!")
        return redirect(request.url)
        # return redirect(url_for('_file', name=filename))
    # return


@app.route('/aws', methods=['POST', 'GET'])
def crispr_aws():
    ec2 = boto3.client('ec2',
                       'us-west-2',
                       aws_access_key_id=request.form.get('pwd'),
                       aws_secret_access_key=request.form.get('key'))
    reservations = ec2.describe_instances(
        Filters=[
            {'Name': 'instance-state-name',
             'Values': [
                 'running',
             ]},
            {'Name': 'image-id',
             'Values': [
                 'ami-07165a1a20b052514',
             ]}])
    if len(reservations.get('Reservations')) == 0:
        print("No existing instance found. Creating a new one!")
        instances = ec2.run_instances(InstanceType="t2.large",
                                      MaxCount=1,
                                      MinCount=1,
                                      ImageId="ami-07165a1a20b052514",
                                      KeyName="snigdha",
                                      IamInstanceProfile={
                                          'Arn': 'arn:aws:iam::423543210473:role/S3fromEC2'
                                      },
                                      BlockDeviceMappings=[
                                          {
                                              'DeviceName': '/dev/sde',
                                              'Ebs': {
                                                  'DeleteOnTermination': True,
                                                  'VolumeSize': 1000,
                                                  'VolumeType': 'standard',
                                                  'Encrypted': False
                                              }
                                          },
                                      ],
                                      TagSpecifications=[{'Tags': [
                                          {
                                              'Key': 'Project',
                                              'Value': 'opencell-sequencing-primer-design'
                                          },
                                          {
                                              'Key': 'Team Leader',
                                              'Value': 'Manuel Leonetti'
                                          },
                                          {
                                              'Key': 'Name',
                                              'Value': 'ML_crispr',
                                          },
                                      ]}])
        public_ip = instances["Instances"][0]["PublicDnsName"]
    else:
        public_ip = reservations.get('Reservations')[0].get('Instances')[0].get('PublicDnsName')

    fd, path = tempfile.mkstemp()
    file = request.files['keyfile']
    os.write(fd, file.read())
    os.system('scp -i %s %s ubuntu@%s://home/ubuntu/' % (path, app.config['CSV'], public_ip))
    os.system('scp -i %s %s ubuntu@%s://home/ubuntu/' % (path, 'crispr_primer.py', public_ip))
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    privkey = paramiko.RSAKey.from_private_key_file(path)
    ssh.connect(public_ip, username='ubuntu', pkey=privkey)
    now = datetime.now()
    timestamp = now.strftime("%m-%d-%Y_%H-%M-%S")
    input_csv = app.config['CSV'].split("/")[1]
    output_name = input_csv.split('_in.csv')[0] + '_out.csv'

    _, stdout, stderr = ssh.exec_command('/home/ubuntu/anaconda/bin/pip install --upgrade pip')
    print(stdout.readlines(), stderr.readlines())
    _, stdout, stderr = ssh.exec_command('/home/ubuntu/anaconda/bin/pip install fastinterval')
    print(stdout.readlines(), stderr.readlines())
    print('nohup /home/ubuntu/anaconda/bin/python crispr_primer.py -f ' + \
                 input_csv + ' -g hg38 -o ' + output_name + ' >> ' + timestamp + '.log &')
    _, stdout, stderr = ssh.exec_command('nohup /home/ubuntu/anaconda/bin/python crispr_primer.py -f ' + \
                 input_csv + ' -g hg38 -o ' + output_name + ' >> ' + timestamp + '.log &')
    # cmd_to_run = '/home/ubuntu/anaconda/bin/pip install --upgrade pip && /home/ubuntu/anaconda/bin/pip install fastinterval && nohup /home/ubuntu/anaconda/bin/python crispr_primer.py -f ' + \
    #              app.config['CSV'] + ' -g hg38 -o ' + output_name + ' >> ' + timestamp + '.log &'
    # _, stdout, stderr = ssh.exec_command(cmd_to_run, timeout=None, get_pty=True)
    # print(stdout.readlines())
    # print(stderr.readlines())
    # exit_status = stdout.channel.recv_exit_status()

    # sftp = ssh.open_sftp()
    # if exit_status == 0:
    #     # try:
    #     #     sftp.stat(output_name)
    #     flash("Primer file generated!")
    # else:
    #     flash("There was an error generating primer file :(")

    return redirect(request.referrer)


# @app.route('/my-link/')
# def my_link():
#   step1.delimit_insertion()

if __name__ == '__main__':
    app.secret_key = '\\\x05\xe4k\x9d\xccq\x99k\x18\xb8\xa58\x96\xcb\x1f'
    app.config['SESSION_TYPE'] = 'filesystem'
    app.debug = True
    app.run(debug=True)
