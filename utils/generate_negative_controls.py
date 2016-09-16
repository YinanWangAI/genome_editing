import numpy as np
import pandas as pd
import regex
import sqlalchemy


def generate_random_sgrna(upstream=20, downstream=0):
    """Generate random sgRNA, the sequence of sgRNA is: upstream NGG downstream

    Args:
        upstream: int, the number of upstream base pairs
        downstream: int, the number of downstream base pairs

    Returns:
        str, a random sgRNA
    """
    bps = np.array(['A', 'T', 'C', 'G'])
    up = np.random.choice(bps, upstream + 1, replace=True)
    down = np.random.choice(bps, downstream, replace=True)
    random_seq = "".join(up) + 'GG' + "".join(down)
    return random_seq


def whole_genome_sgrna(
        upstream=20, downstream=0, only_upstream=True,
        uri='postgresql://yinan:123456@localhost/genome_editing'):
    """All sgRNAs in Grch38 genome

    Args:
        upstream: int, the number of upstream base pairs
        downstream: int, the number of downstream base pairs
        only_upstream: boolean, whether only return upstream sequence
        uri: sqlalchemy URI

    Returns:
        [gg_seq, cc_seq],
        gg_seq: all sgRNAs with sequence: upstream NGG downstream
        cc_seq: all sgRNAs with sequence: downstream CCN upstream
    """
    # whole genome chromosomes
    chroms = ['chr' + str(x) for x in range(1, 23)]
    chroms.append('chrX')
    chroms.append('chrY')
    chroms.append('chrM')

    # database engine
    engine = sqlalchemy.create_engine(uri)

    # get all sgRNAs
    gg_seq = []
    cc_seq = []
    for chrom in chroms:
        print(chrom)
        table_name = 'grch38_' + chrom
        chr_seq = pd.read_sql(table_name, engine).iloc[0, 0].upper()
        gg_match = regex.finditer('\w{%d}\wGG\w{%d}' % (upstream, downstream),
                                  chr_seq, overlapped=True)
        cc_match = regex.finditer('\w{%d}CC\w\w{%d}' % (downstream, upstream),
                                  chr_seq, overlapped=True)
        if only_upstream:
            for sub_seq in gg_match:
                gg_seq.append(sub_seq.group()[:upstream])
            for sub_seq in cc_match:
                cc_seq.append(sub_seq.group()[-upstream:])
        else:
            for sub_seq in gg_match:
                gg_seq.append(sub_seq.group())
            for sub_seq in cc_match:
                cc_seq.append(sub_seq.group())
    print('Convert to numpy array')
    gg_seq = np.asarray(gg_seq)
    cc_seq = np.asarray(cc_seq)
    return [gg_seq, cc_seq]


def cut_seq(ref_seqs, cut_num=3, length=None, single_seq=False):
    """Cut sequence into non-overlap fragments

    Args:
        ref_seqs: str or numpy array, the sequence to be cut
        cut_num: the number of fragments
        length: the length of sequence
        single_seq: whether ref_seqs is a single sequence or an array

    Returns:
        cut_seqs, a dict contains the sequences for each fragment
    """
    if single_seq:
        ref_seqs = [ref_seqs]

    if length is None:
        seq_len = len(ref_seqs[0])
    else:
        seq_len = length

    fragment_size = int(np.ceil(seq_len / cut_num))
    cut_seqs = {}
    for i in range(cut_num):
        start = i * fragment_size
        end = min((i + 1) * fragment_size, seq_len)
        cut_seqs[i] = [x[start:end] for x in ref_seqs]
    return cut_seqs
