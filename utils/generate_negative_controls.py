from Bio.Seq import Seq
import numpy as np
import multiprocessing as mp
import pandas as pd
import pickle
import regex
import sqlalchemy


def generate_negative_control_mp(num, file_path, whole_genome_seq=None,
                                 length=20, min_mismatch=2, seed=1015,
                                 num_cpu=4):
    p = mp.Pool(num_cpu)
    for i in range(num_cpu):
        p.apply_async(generate_negative_control,
                      args=(num, whole_genome_seq, length, min_mismatch, seed,
                            file_path))
    p.close()
    p.join()


def generate_negative_control(num, whole_genome_seq=None, length=20,
                              min_mismatch=2, seed=None, file_path=None):
    # set seed
    np.random.seed(seed)

    # generate whole genome sgRNAs
    print('Preparing whole genome sgRNA references')
    if whole_genome_seq is None:
        whole_genome_seq = whole_genome_sgrna(upstream=length,
                                              only_upstream=True)
        gg_seq_cut = whole_genome_seq[0]
        cc_seq_cut = whole_genome_seq[1]
    else:
        gg_seq_cut = whole_genome_seq[0]
        cc_seq_cut = whole_genome_seq[1]
    del whole_genome_seq

    # open file which stores results
    if file_path is not None:
        file_handler = open(file_path, 'a')

    # generate negative controls
    neg_controls = []
    while len(neg_controls) < num:
        print('Generating sgRNA')
        random_seq = generate_random_sgrna(upstream=length, downstream=0)[:-3]

        # remove sgRNAs containing TTTT
        if random_seq.find('TTTT') != -1:
            continue
        random_seq_cut = list(random_seq)

        # compare with *NGG*
        print('Comparing with *NGG*')
        if is_negative_control(gg_seq_cut, random_seq_cut, min_mismatch):
            # compare with *CCN*
            print('Comparing with *CCN*')
            rc_seq = str(Seq(random_seq).reverse_complement())
            random_rc_seq_cut = list(rc_seq)
            if is_negative_control(cc_seq_cut, random_rc_seq_cut, min_mismatch):
                neg_controls.append(random_seq)
                file_handler.write(random_seq + '\n')
    file_handler.close()
    return list(neg_controls)


def whole_genome_sgrna(
        upstream=20, downstream=0, only_upstream=True,
        uri='postgresql://wangyn@localhost/genome_editing'):
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
                gg_seq.append(list(sub_seq.group()[:upstream]))
            for sub_seq in cc_match:
                cc_seq.append(list(sub_seq.group()[-upstream:]))
        else:
            for sub_seq in gg_match:
                gg_seq.append(list(sub_seq.group()))
            for sub_seq in cc_match:
                cc_seq.append(list(sub_seq.group()))
    print('Convert to numpy array')
    gg_seq = np.asarray(gg_seq)
    cc_seq = np.asarray(cc_seq)
    return [gg_seq, cc_seq]


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


def is_negative_control(cut_seq_1, cut_seq_2, min_mismatch):
    num_mismatch = np.min(np.sum(cut_seq_1 != cut_seq_2, axis=1))
    if num_mismatch > min_mismatch:
        return True
    else:
        return False


if __name__ == '__main__':
    with open('./whole_genome_sgrna_mat.pkl', 'rb') as f:
        whole_genome_sgrna_mat = pickle.load(f)

