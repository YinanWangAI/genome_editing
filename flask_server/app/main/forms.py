from flask_wtf import Form
from wtforms import StringField, SubmitField, IntegerField, SelectField, \
    TextAreaField
from wtforms.validators import DataRequired


class DesignSgrnaForm(Form):
    design_input = TextAreaField('Gene Symbols/Refseq IDs/Sequences',
                                 validators=[DataRequired()],
                                 render_kw={'rows': 18})
    pams = SelectField('PAM', choices=[('NGG', 'NGG'), ('NAG', 'NAG')],
                       validators=[DataRequired()])
    ref_genome = SelectField('Reference Genome',
                             choices=[('hg19', 'hg19'), ('hg38', 'hg38'),
                                      ('mm10', 'mm10')],
                             validators=[DataRequired()])
    # upstream_len = IntegerField('Upstream Length',
    #                             render_kw={'value': 4},
    #                             validators=[DataRequired()])
    # downstream_len = IntegerField('Downstream Length',
    #                               render_kw={'value': 3},
    #                               validators=[DataRequired()])
    sgrna_len = IntegerField('sgRNA Length',
                             render_kw={'value': 20},
                             validators=[DataRequired()])
    flank_len = IntegerField('Flank Length',
                             render_kw={'value': 30},
                             validators=[DataRequired()])
    # filter_tttt = SelectField('Filter TTTT?',
    #                           choices=[('Yes', 'Yes'), ('No', 'No')])
    only_target_aa = SelectField('Only Target Amino Acid?',
                                 choices=[('Yes', 'Yes'), ('No', 'No')],
                                 validators=[DataRequired()])
    input_type = SelectField('Input Type',
                             choices=[('Gene Symbol', 'Gene Symbol'),
                                      ('Refseq ID', 'Refseq ID'),
                                      ('Sequence', 'Sequence')],
                             validators=[DataRequired()])
    submit = SubmitField('Submit')


class ScoreSgrnaForm(Form):
    seqs = TextAreaField('sgRNA Sequences', validators=[DataRequired()],
                         render_kw={'rows': 5})
    score_algo = SelectField('Score Algorithm',
                             choices=[('rs2', 'rs2')])
    submit = SubmitField('Submit')
