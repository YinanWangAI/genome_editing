import os
from genome_editing.utils.alignment import bowtie_alignment

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
    seq_list = [seq]
    for char in pam:
        seq_list = [seq + x for seq in seq_list for x in convert_dict[char]]
    return seq_list


def have_off_targets(seq, pam, upstream_len=16, num_mismatch=0,
                        bowtie_index=HG38_BOWTIE_INDEX_PATH):
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


if __name__ == '__main__':
    seq = 'ATAAGCGATTTTCCACTGGG'
    pam = 'NGG'
    upstream_len = 16
    num_mismatch = 0
    have_off_targets(seq, pam, upstream_len=16, num_mismatch=0,
                     bowtie_index=HG38_BOWTIE_INDEX_PATH)
