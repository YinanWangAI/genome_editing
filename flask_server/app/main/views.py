import sys
sys.path.append('/Users/yinan/PycharmProjects/')
import genome_editing.design_sgRNA.design as dsr
from flask import Flask, render_template, redirect, url_for, send_from_directory
from . import main
from .forms import DesignSgrnaForm
from .. import db


JOB_ID = 0

@main.route('/')
def main_page():
    return render_template("main_page.html")


@main.route('/design_sgrna/', methods=['GET', 'POST'])
def design_sgrna():
    global JOB_ID
    form = DesignSgrnaForm()
    if form.validate_on_submit():
        JOB_ID += 1
        gene_ids = form.gene_ids.data.split('\n')
        upstream_len = form.upstream_len.data
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


@main.route('/design_sgrna/<job_id>/')
def design_sgrna_download_link(job_id):
    result_dir = '/Users/yinan/PycharmProjects/genome_editing/flask_server/results'
    file_name = 'design_sgrna_output_' + str(job_id) + '.csv'
    return send_from_directory(result_dir, file_name, as_attachment=True)


@main.route('/under_development')
def under_development():
    return render_template('under_development.html')