import numpy as np
import os
import pandas as pd
import subprocess
import tempfile
from genome_editing.utils.alignment import bowtie_alignment
from genome_editing.utils.utilities import reverse_complement

"""
1. map the sequence without PAM to genome
2. if no target except the target gene, pass
3. if target: extract the position then the sequence, check whether the
  off-target is on the upstream of PAM
"""

HG38_BOWTIE_INDEX_PATH = os.getenv('HG38_BOWTIE_INDEX_PATH')


def extend_seq(seq, pam):
    """Concat sequence and PAM"""
    convert_dict = {'A': 'A', 'T': 'T', 'C': 'C', 'G': 'G', 'N': 'ATCG'}
    if type(seq) == str:
        seq_list = [seq]
    else:
        seq_list = seq
    for char in pam:
        seq_list = [seq + x for seq in seq_list for x in convert_dict[char]]
    return seq_list


def have_off_targets(seq, pam, upstream_len=16, num_mismatch=1,
                     bowtie_index=HG38_BOWTIE_INDEX_PATH):
    """num_mismatch = 1 because there's a N in PAM"""
    trim_seq = seq[-upstream_len:].upper()
    trim_seq_w_pam = extend_seq(trim_seq, pam)
    off_target = False
    for input_seq in trim_seq_w_pam:
        alignment_out = bowtie_alignment(seq=input_seq,
                                         bowtie_index_path=bowtie_index,
                                         num_mismatch=num_mismatch,
                                         seed=upstream_len)
        if alignment_out.shape[0] > 1:
            off_target = True
            break
    return off_target


def sgrna_off_targets(seq, pam='NGG', seed=20, num_mismatch=1,
                      bowtie_index=HG38_BOWTIE_INDEX_PATH):
    off_target = False
    sub_seq = seq[-seed:] + pam
    alignment_out = sgrna_alignment(seq=sub_seq,
                                     bowtie_index_path=bowtie_index,
                                     num_mismatch=num_mismatch)
    if alignment_out.shape[0] > 1:
        off_target = True
    return off_target


def sgrna_alignment(seq, bowtie_index_path=HG38_BOWTIE_INDEX_PATH,
                    num_mismatch=1):
    """num_mismatch = 1 because there's a N in PAM"""
    temp_file = tempfile.NamedTemporaryFile()
    alignment_out_path = temp_file.name

    cmd = 'bowtie -a -p 4 -v {} {} -c {} -S {}'.format(
        num_mismatch, bowtie_index_path, seq, alignment_out_path)

    run_code = subprocess.call(cmd, shell=True)
    if run_code == 0:
        alignment_out = pd.read_table(alignment_out_path,
                                      comment='@', header=None)
        return alignment_out
    else:
        print('Error in Bowtie')
        return 1


def sgrna_off_targets_batch(seqs, pam='NGG', seed=20, num_mismatch=1,
                            bowtie_index=HG38_BOWTIE_INDEX_PATH):
    sub_seqs = [seq[-seed:] + pam for seq in seqs]
    alignment_out = sgrna_alignment_batch(seqs=sub_seqs,
                                          bowtie_index_path=bowtie_index,
                                          num_mismatch=num_mismatch)
    alignment_out_seqs = alignment_out.iloc[:, 4]
    off_targets = []
    for seq in sub_seqs:
        target_num = np.sum(alignment_out_seqs == seq) + \
                     np.sum(alignment_out_seqs == reverse_complement(seq))
        if target_num > 1:
            off_targets.append(True)
        else:
            off_targets.append(False)
    return off_targets


def sgrna_alignment_batch(seqs, bowtie_index_path=HG38_BOWTIE_INDEX_PATH,
                          num_mismatch=1):
    """num_mismatch = 1 because there's a N in PAM"""
    temp_file = tempfile.NamedTemporaryFile()
    alignment_out_path = temp_file.name

    cmd = 'bowtie -a -p 4 -v {} {} -c {} {}'.format(
        num_mismatch, bowtie_index_path, ','.join(seqs), alignment_out_path)

    run_code = subprocess.call(cmd, shell=True)
    if run_code == 0:
        alignment_out = pd.read_table(alignment_out_path,
                                      comment='@', header=None)
        return alignment_out
    else:
        print('Error in Bowtie')
        return 1


if __name__ == '__main__':
    seq = 'ATAAGCGATTTTCCACTGGG'
    pam = 'NGG'
    upstream_len = 16
    num_mismatch = 0
    have_off_targets(seq, pam, upstream_len=16, num_mismatch=0,
                     bowtie_index=HG38_BOWTIE_INDEX_PATH)
