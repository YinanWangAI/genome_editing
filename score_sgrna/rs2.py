"""Rule set 2 algorithm from:

Optimized sgRNA design to maximize activity and minimize off-target effects of
CRISPR-Cas9. Nature Biotechnology, 1–12. http://doi.org/10.1038/nbt.3437
"""
import os
import numpy as np
import pandas as pd
import subprocess

PYTHON2 = os.getenv('ANACONDA_PYTHON2')
RS2 = os.getenv('RS2_CALCULATOR')


def compute_rs2(seq, aa_cut=None, per_peptide=None):
    if aa_cut is None:
        aa_cut = -1
    if per_peptide is None:
        per_peptide = -1
    cmd = PYTHON2 + ' ' + RS2 + ' --seq {} --aa-cut {} --per-peptide {}'.format(
        seq, aa_cut, per_peptide)
    rs2_score = subprocess.check_output(cmd, shell=True).decode('utf-8')
    return float(rs2_score.strip().split(' ')[-1])


def compute_rs2_batch(seqs):
    if type(seqs) is not list:
        seqs = [seqs]
    flag = True
    seq_num = 0
    for seq in seqs:
        if (seq == '') or (len(seq) < 30):
            continue
        seq = seq[:30]
        seq_score = compute_rs2(seq)
        if flag:
            seq_num += 1
            out = np.asarray([seq, seq_score])
            flag = False
        else:
            seq_num += 1
            temp = np.asarray([seq, seq_score])
            out = np.vstack((out, temp))
    if seq_num < 2:
        out = pd.DataFrame(out).transpose()
    else:
        out = pd.DataFrame(out)
    out.columns = ['seq', 'rs2_score']
    return out


if __name__ == '__main__':
    print(compute_rs2('CAGAAAAAAAAACACTGCAACAAGAGGGTA'))
