from flask import Flask, render_template, redirect, url_for, send_from_directory
from flask_bootstrap import Bootstrap
from flask_migrate import Migrate, MigrateCommand
from flask_moment import Moment
from flask_script import Manager
from flask_sqlalchemy import SQLAlchemy
from flask_wtf import Form
from wtforms import StringField, SubmitField, IntegerField, SelectField, \
    TextAreaField
from wtforms.validators import data_required

import sys
sys.path.append('/Users/yinan/PycharmProjects/')
import genome_editing.design_sgRNA.design as dsr


JOB_ID = 0


app = Flask(__name__)
app.config['SECRET_KEY'] = 'secret key'
app.config['SQLALCHEMY_DATABASE_URI'] = \
    'postgresql://yinan:123456@localhost/genome_editing'
app.config['SQLALCHEMY_COMMIT_ON_TEARDOWN'] = True
bootstrap = Bootstrap(app)
manager = Manager(app)
moment = Moment(app)
db = SQLAlchemy(app)
migrate = Migrate(app, db)
manager.add_command('db', MigrateCommand)


@app.route('/')
def main_page():
    return render_template("main_page.html")


class DesignSgrnaForm(Form):
    gene_ids = TextAreaField('Entrez IDs', validators=[data_required()],
                             render_kw={'rows':3})
    pams = StringField('PAM', validators=[data_required()],
                       render_kw={'value':'NGG'})
    upstream_len = IntegerField('Upstream Length', validators=[data_required()],
                                render_kw={'value': 4})
    downstream_len = IntegerField('Downstream Length',
                                  validators=[data_required()],
                                  render_kw={'value': 3})
    sgrna_len = IntegerField('sgRNA Length', validators=[data_required()],
                             render_kw={'value': 20})
    flank_len = IntegerField('Flank Length', validators=[data_required()],
                             render_kw={'value': 30})
    filter_tttt = SelectField('Filter TTTT?',
                              choices=[('Yes', 'Yes'), ('No', 'No')])
    submit = SubmitField('Submit')


# class DesignResults(db.Model):
#     __tablename__ = 'design_reuslts'
#     gene_id = db.Column(db.Integer, index=True)
#     sgrna_seq = db.Column(db.String, primary_key=True)
#
#     def __repr__(self):
#         return 'gene: {}'.format(self.gene_id)
#
#
# class GeneInfo(db.Model):
#     __tablename__ = 'gene_info'
#     gene_id = db.Column(db.Integer, db.ForeignKey('design_results.gene_id'))
#     gene_symbol = db.Column(db.String)
#
#     def __repr__(self):
#         return 'gene: {}'.format(self.gene_symbol)


@app.route('/design_sgrna/', methods=['GET', 'POST'])
def design_sgrna():
    global JOB_ID
    form = DesignSgrnaForm()
    if form.validate_on_submit():
        JOB_ID += 1
        gene_ids = form.gene_ids.data.split('\n')
        upstream_len =form.upstream_len.data
        downstream_len = form.downstream_len.data
        flank_len = form.flank_len.data
        sgrna_len = form.sgrna_len.data
        pam = form.pams.data
        if form.filter_tttt.data == 'Yes':
            filter_tttt = True
        else:
            filter_tttt = False

        flag = True
        for gene_id in gene_ids:
            gene_id = int(gene_id)
            designer = dsr.Designer(entrez_id=gene_id,
                                    sgrna_upstream=upstream_len,
                                    sgrna_downstream=downstream_len,
                                    sgrna_length=sgrna_len,
                                    flank=flank_len,
                                    overlapped=True,
                                    filter_tttt=filter_tttt)
            designer.get_sgrnas(pams=[pam])
            if flag:
                sgrna_designer_out = designer.output()
                flag = False
            else:
                sgrna_designer_out = sgrna_designer_out.append(
                    designer.output())
        output_path = './results/design_sgrna_output_' + str(JOB_ID) + '.csv'
        sgrna_designer_out.to_csv(output_path, index=None)
        return render_template('design_sgrna_output.html', job_id=JOB_ID)
    else:
        return render_template('design_sgrna_input.html', form=form)


@app.route('/design_sgrna/<job_id>/')
def design_sgrna_download_link(job_id):
    result_dir = '/Users/yinan/PycharmProjects/genome_editing/flask_server/results'
    file_name = 'design_sgrna_output_' + str(job_id) + '.csv'
    return send_from_directory(result_dir, file_name, as_attachment=True)


@app.route('/under_development')
def under_development():
    return render_template('under_development.html')


@app.errorhandler(404)
def page_not_found(e):
    return render_template('404.html'), 404


@app.errorhandler(500)
def internal_server_error(e):
    return render_template('500.html'), 500


if __name__ == "__main__":
    manager.run()
