from flask_wtf import Form
from wtforms import StringField, SubmitField, IntegerField, SelectField, \
    TextAreaField
from wtforms.validators import data_required


class DesignSgrnaForm(Form):
    gene_ids = TextAreaField('Entrez IDs', validators=[data_required()],
                             render_kw={'rows':3})
    pams = StringField('PAM', validators=[data_required()],
                       render_kw={'value':'NGG'})
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
    submit = SubmitField('Submit')
