import numpy as np
import pandas as pd
import regex


def get_reads_info(fq_1, fq_2=None, pattern='ACCG.*GTTTA', quick_merge=True):
    with open(fq_1) as f:
        sgrnas_1 = []
        for i, line in enumerate(f):
            if i % 4 == 0:
                record = line.split(' ')[0]
            if i % 4 == 1:
                seq = regex.findall(pattern, line)
                if len(seq) == 0:
                    sgrnas_1.append([record, np.nan])
                else:
                    sgrnas_1.append([record, seq[0]])
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
                    seq = regex.findall(pattern, line)
                    if len(seq) == 0:
                        sgrnas_2.append([record, np.nan])
                    else:
                        sgrnas_2.append([record, seq[0]])
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
