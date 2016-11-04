from flask_wtf import Form
from wtforms import StringField, SubmitField, IntegerField, SelectField, \
    TextAreaField
from wtforms.validators import DataRequired


class DesignSgrnaBase(Form):
    ref_genome = SelectField('Reference Genome',
                             choices=[('hg19', 'hg19'), ('hg38', 'hg38'),
                                      ('mm10', 'mm10')],
                             validators=[DataRequired()])
    input_type = SelectField('Input Type',
                             choices=[('Gene Symbol', 'Gene Symbol'),
                                      ('Refseq ID', 'Refseq ID'),
                                      ('Sequence', 'Sequence')],
                             validators=[DataRequired()])
    submit = SubmitField('Submit')


class DesignBatchSgrnaForm(DesignSgrnaBase):
    design_input = TextAreaField('Gene Symbols/Refseq IDs',
                                 validators=[DataRequired()],
                                 render_kw={'rows': 10})
    pams = SelectField('PAM', choices=[('NGG', 'NGG'), ('NAG', 'NAG')],
                       validators=[DataRequired()])


class DesignSingleSgrnaForm(DesignSgrnaBase):
    design_input = TextAreaField('Gene Symbol/Refseq ID/Sequence',
                                 render_kw={'rows': 10},
                                 validators=[DataRequired()])
    pam_seq = StringField('PAM Sequence', default='NGG',
                          validators=[DataRequired()])
    upstream_len = IntegerField('Upstream Length', default=4,
                                render_kw={'min': '0', 'step': '1',
                                           'type': 'number'})
    downstream_len = IntegerField('Downstream Length', default=3,
                                  render_kw={'min': '0', 'step': '1',
                                             'type': 'number'})
    sgrna_len = IntegerField('sgRNA Length', default=20,
                             render_kw={'min': '1', 'step': '1',
                                        'type': 'number'},
                             validators=[DataRequired()])
    flank_len = IntegerField('Flank Length', default=30,
                             render_kw={'min': '0', 'step': '1',
                                        'type': 'number'})
    filter_tttt = SelectField('Filter TTTT?',
                              choices=[('Yes', 'Yes'), ('No', 'No')])


class ScoreSgrnaForm(Form):
    seqs = TextAreaField('sgRNA Sequences', validators=[DataRequired()],
                         render_kw={'rows': 5})
    score_algo = SelectField('Score Algorithm',
                             choices=[('rs2', 'rs2')])
    submit = SubmitField('Submit')
