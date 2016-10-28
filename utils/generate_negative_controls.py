import numpy as np
from genome_editing.utils import alignment


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


def generate_neg_control(num, length, pam='', off_targets=2,
                         seed=None, file_path=None):
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
        random_seq += pam
        alignment_info = alignment.bowtie_alignment(seq=random_seq,
                                                    report_all=False)
        if alignment_info.iloc[:, 1].values == 4:
            if random_seq[:length] not in neg_controls:
                if f:
                    f.write(random_seq[:length] + '\n')
                neg_controls.append(random_seq[:length])
                if len(neg_controls) % 100 == 0:
                    print('Generate {} negative controls'.format(
                        len(neg_controls)))
    if f:
        f.close()
    return neg_controls


if __name__ == '__main__':
    print(generate_neg_control(num=10, length=20, pam=''))
