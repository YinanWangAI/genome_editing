"""Compute SCC score from GR2015"""

import os
import pandas as pd
import subprocess
import tempfile


def compute_scc(seqs,
                scc_path='/Users/yinan/PycharmProjects/genome_editing/score_sgrna/SSC0.1/bin/SSC',
                mat_path='/Users/yinan/PycharmProjects/genome_editing/score_sgrna/SSC0.1/matrix/human_mouse_CRISPR_KO_30bp.matrix'):
    """Compute SCC score

    Args:
        seqs: list, input sequence, 20mer + PAM + 7mer
        scc_path: the path of SCC bin
        mat_path: the path of SCC score matrix

    Returns:
        DataFrame, seqs and SSC score
    """
    temp_input = tempfile.NamedTemporaryFile(mode='a', delete=False)
    temp_output = tempfile.NamedTemporaryFile(mode='a', delete=False)
    for seq in seqs:
        temp_input.write(seq + '\n')
    temp_input.close()
    cmd = '{} -l 30 -m {} -i {} -o {}'.format(scc_path, mat_path,
                                              temp_input.name, temp_output.name)
    temp_output.close()
    subprocess.run(cmd, shell=True)
    scc = pd.read_table(temp_output.name, header=None)
    scc.columns = ['seq_with_context', 'scc_score']
    os.remove(temp_input.name)
    os.remove(temp_output.name)
    return scc


def compute_scc_crispr_ia(seqs, spacer_len,
                scc_path='/Users/yinan/PycharmProjects/genome_editing/score_sgrna/SSC0.1/bin/SSC',
                mat_path_prefix='/Users/yinan/PycharmProjects/genome_editing/score_sgrna/SSC0.1/matrix/'):
    temp_input = tempfile.NamedTemporaryFile(mode='a', delete=False)
    temp_output = tempfile.NamedTemporaryFile(mode='a', delete=False)
    for seq in seqs:
        temp_input.write(seq + '\n')
    temp_input.close()
    mat_path = mat_path_prefix + 'human_CRISPRi_{}bp.matrix'.format(spacer_len)
    cmd = '{} -l {} -m {} -i {} -o {}'.format(scc_path, spacer_len, mat_path,
                                              temp_input.name, temp_output.name)
    temp_output.close()
    subprocess.run(cmd, shell=True)
    scc = pd.read_table(temp_output.name, header=None)
    scc.columns = ['spacer_seq', 'scc_score']
    os.remove(temp_input.name)
    os.remove(temp_output.name)
    return scc
