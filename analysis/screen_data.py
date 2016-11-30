import numpy as np
import pandas as pd
import regex


def get_reads_info(fq_1, fq_2=None, pattern='ACCG.{20}GTTTA', quick_merge=True):
    pattern = regex.compile(pattern)
    with open(fq_1) as f:
        sgrnas_1 = []
        for i, line in enumerate(f):
            if i % 4 == 0:
                record = line.split(' ')[0]
            if i % 4 == 1:
                seq = pattern.search(line)
                if seq is None:
                    sgrnas_1.append([record, np.nan])
                else:
                    sgrnas_1.append([record, seq.captures()[0]])
            if i % 10000000 == 0:
                print(i)
        sgrnas_df = pd.DataFrame(sgrnas_1)
        sgrnas_df.columns = ['reads_id', 'seq_1']

    if fq_2 is not None:
        with open(fq_2) as f:
            sgrnas_2 = []
            for i, line in enumerate(f):
                if i % 4 == 0:
                    record = line.split(' ')[0]
                if i % 4 == 1:
                    seq = pattern.search(line)
                    if seq is None:
                        sgrnas_2.append([record, np.nan])
                    else:
                        sgrnas_2.append([record, seq.captures()[0]])
                if i % 10000000 == 0:
                    print(i)
            sgrnas_df_2 = pd.DataFrame(sgrnas_2)
            sgrnas_df_2.columns = ['reads_id', 'seq_2']
        if quick_merge:
            print('quick merge...')
            rand_index = np.random.choice(sgrnas_df.shape[0], 100,
                                          replace=False)
            reads_id1 = sgrnas_df.iloc[rand_index, :].reads_id.values
            reads_id2 = sgrnas_df_2.iloc[rand_index, :].reads_id.values
            if np.all(reads_id1 == reads_id2):
                sgrnas_df.loc[:, 'seq_2'] = sgrnas_df_2.loc[:, 'seq_2'].values
            else:
                print('Fail to quick merge, transfer to standard merge')
                sgrnas_df = pd.merge(sgrnas_df, sgrnas_df_2, on='reads_id')
        else:
            sgrnas_df = pd.merge(sgrnas_df, sgrnas_df_2, on='reads_id')
    return sgrnas_df


def decode_summary(df):
    reads_1 = df.loc[:, ['reads_id', 'seq_1']].dropna().reads_id.values
    reads_2 = df.loc[:, ['reads_id', 'seq_2']].dropna().reads_id.values
    intersect_reads = set(reads_1) & set(reads_2)
    union_reads = set(reads_1) | set(reads_2)

    double_map_ratio = len(intersect_reads) / len(union_reads)
    map_ratio = len(union_reads) / df.shape[0]
    return map_ratio, double_map_ratio


def count_sgrna(df, start=0, end=-1):
    """统计reads数

    Args:
        df: Output of get_reads_info
        start: start index of seq
        end: end index of seq
    Returns:

    """
    df = df[pd.isnull(df.seq_1) | pd.isnull(df.seq_2.values)]
    df_seqs = list(df.seq_1.dropna().values) + list(df.seq_2.dropna().values)
    count_dict = {}
    for seq in df_seqs:
        if seq not in count_dict:
            count_dict[seq] = 1
        else:
            count_dict[seq] += 1
    reads_count = pd.DataFrame.from_dict(count_dict, orient='index')
    reads_count.loc[:, 'sgrna_seq'] = [x[start:end] for x in
                                       list(reads_count.index)]
    reads_count.index = range(reads_count.shape[0])
    reads_count.columns = ['counts', 'sgrna_seq']
    return reads_count
