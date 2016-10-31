import sys
sys.path.append('/Users/yinan/PycharmProjects/')
import genome_editing.design_sgRNA.design as dsr
from genome_editing.score_sgrna.rs2 import compute_rs2_batch
from flask import Flask, render_template, redirect, url_for, send_from_directory
from . import main
from .forms import DesignSgrnaForm, ScoreSgrnaForm
from .. import db
from ..models import SgrnaDesign
from genome_editing.utils.utilities import model_to_df, coordinate_sgrna


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

        input_type = form.input_type.data
        design_inputs = form.design_input.data
        # upstream_len = form.upstream_len.data
        # downstream_len = form.downstream_len.data
        flank_len = form.flank_len.data
        sgrna_len = form.sgrna_len.data
        pam = form.pams.data

        # TODO: support different up-, down- stream and sgRNA length
        if input_type == 'Gene Symbol':
            gene_symbols = design_inputs.split('\n')
            gene_symbols = [x.strip() for x in gene_symbols]
            print(gene_symbols)
            query_out = SgrnaDesign.query.filter(
                SgrnaDesign.gene_symbol.in_(gene_symbols)).all()
            sgrna_designer_out = model_to_df(query_out)
            # sgrna_designer_out = coordinate_sgrna(
            #     sgrna_designer_out, upstream_len, downstream_len, sgrna_len)
        elif input_type == 'Refseq ID':
            refseq_ids = design_inputs.split('\n')
            refseq_ids = [x.strip() for x in refseq_ids]
            print(refseq_ids)
            query_out = SgrnaDesign.query.filter(
                SgrnaDesign.refseq_id.in_(refseq_ids)).all()
            sgrna_designer_out = model_to_df(query_out)
            # sgrna_designer_out = coordinate_sgrna(
            #     sgrna_designer_out, upstream_len, downstream_len, sgrna_len)
        else:
            seq = design_inputs.split('\n')
            assert len(seq) == 1, "Too many sequences"
            seq = seq[0]
            sgrna_design = dsr.SeqDesigner(seq=seq,
                                           sgrna_upstream=upstream_len,
                                           sgrna_downstream=downstream_len,
                                           sgrna_length=sgrna_len,
                                           flank=flank_len)
            sgrna_design.get_sgrnas(pams=[pam])
            sgrna_designer_out = sgrna_design.output()

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


@main.route('/score_sgrna/', methods=['GET', 'POST'])
def score_sgrna():
    global JOB_ID
    form = ScoreSgrnaForm()

    if form.validate_on_submit():
        JOB_ID += 1
        seqs = form.seqs.data.split('\n')
        print(seqs)
        score_sgrna_out = compute_rs2_batch(seqs)
        output_path = './results/score_sgrna_output_' + str(JOB_ID) + '.csv'
        score_sgrna_out.to_csv(output_path, index=None)
        return render_template('score_sgrna_output.html', job_id=JOB_ID)
    else:
        return render_template('score_sgrna_input.html', form=form)


@main.route('/score_sgrna/<job_id>/')
def score_sgrna_download_link(job_id):
    result_dir = '/Users/yinan/PycharmProjects/genome_editing/flask_server/results'
    file_name = 'score_sgrna_output_' + str(job_id) + '.csv'
    return send_from_directory(result_dir, file_name, as_attachment=True)


@main.route('/under_development')
def under_development():
    return render_template('under_development.html')
