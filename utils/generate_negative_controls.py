import os
import numpy as np
from genome_editing.utils import alignment


HG38_BOWTIE_INDEX_PATH = os.getenv('HG38_BOWTIE_INDEX_PATH')


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


def generate_neg_control(num, length, seed_len, num_mismatch=2,
                         seed=None, file_path=None,
                         bowtie_index_path=HG38_BOWTIE_INDEX_PATH):
    """Generate negative controls. No totally match of the seed_len upstream of
    PAM on the whole genome.

    Args:
        num:
        length:
        seed_len:
        num_mismatch:
        seed:
        file_path:
        bowtie_index_path:

    Returns:

    """
    # set seed
    np.random.seed(seed)

    # open file
    if file_path is not None:
        f = open(file_path, 'a')
    else:
        f = None

    # generate negative controls
    neg_controls = []
    while len(neg_controls) < num:
        random_seq = generate_random_sgrna(upstream=length,
                                           downstream=0)[:length]

        # remove sgRNAs containing TTTT
        if random_seq.find('TTTT') != -1:
            continue

        # alignment
        input_seq= random_seq[-seed_len:]
        alignment_info = alignment.bowtie_alignment(
            seq=input_seq, report_all=False, num_mismatch=num_mismatch,
            bowtie_index_path=bowtie_index_path, seed=seed_len)
        if alignment_info.iloc[:, 1].values == 4:
            if random_seq not in neg_controls:
                if f:
                    f.write(random_seq + '\n')
                neg_controls.append(random_seq)
                if len(neg_controls) % 100 == 0:
                    print('Generate {} negative controls'.format(
                        len(neg_controls)))
    if f:
        f.close()
    return neg_controls
