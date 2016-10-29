from flask_wtf import Form
from wtforms import StringField, SubmitField, IntegerField, SelectField, \
    TextAreaField
from wtforms.validators import data_required


class DesignSgrnaForm(Form):
    design_input = TextAreaField('Gene Symbols/Refseq IDs/Sequences',
                                 validators=[data_required()],
                                 render_kw={'rows': 18})
    pams = SelectField('PAM', choices=[('NGG', 'NGG'), ('NAG', 'NAG')])
    upstream_len = IntegerField('Upstream Length',
                                render_kw={'value': 4})
    downstream_len = IntegerField('Downstream Length',
                                  render_kw={'value': 3})
    sgrna_len = IntegerField('sgRNA Length',
                             render_kw={'value': 20})
    flank_len = IntegerField('Flank Length',
                             render_kw={'value': 30})
    filter_tttt = SelectField('Filter TTTT?',
                              choices=[('Yes', 'Yes'), ('No', 'No')])
    only_target_aa = SelectField('Only Target Amino Acid?',
                                 choices=[('Yes', 'Yes'), ('No', 'No')])
    input_type = SelectField('Input Type',
                             choices=[('Gene Symbol', 'Gene Symbol'),
                                      ('Refseq ID', 'Refseq ID'),
                                      ('Sequence', 'Sequence')])
    submit = SubmitField('Submit')


class ScoreSgrnaForm(Form):
    seqs = TextAreaField('sgRNA Sequences', validators=[data_required()],
                         render_kw={'rows':5})
    score_algo = SelectField('Score Algorithm',
                             choices=[('rs2', 'rs2')])
    submit = SubmitField('Submit')
