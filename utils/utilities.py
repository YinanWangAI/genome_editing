from Bio.Seq import Seq
import numpy as np
import pandas as pd
import sqlalchemy
import os

GENOME_EDITING_URI = os.environ.get('GENOME_EDITING_URI')
UPSTREAM = 4
DOWNSTREAM = 3
SGRNA_LEN = 20


def gene_symbol_to_refseq(
        genes,
        table_name='igenome_ucsc_hg19_refgene',
        uri=GENOME_EDITING_URI):
    """Convert gene symbol to Refseq ID

    Args:
        genes: gene symbol
        table_name: table storing gene info
        uri: sqla URI

    Returns:
        Refseq IDs
    """
    engine = sqlalchemy.create_engine(uri)
    gene_refseq = {}
    for gene in genes:
        query = "SELECT * FROM {} WHERE name2='{}'".format(table_name, gene)
        gene_info = pd.read_sql_query(query, engine).drop_duplicates()
        gene_refseq[gene] = gene_info.name.values
    return gene_refseq


def model_to_df(query_out):
    flag = True
    for query_row in query_out:
        df_row = np.array([query_row.gene_symbol, query_row.refseq_id,
                           query_row.exon_id, query_row.chrom, query_row.start,
                           query_row.end, query_row.raw_sequence,
                           query_row.pam_type,
                           query_row.cutting_site_type, query_row.cutting_site,
                           query_row.sgrna_seq, query_row.sgrna_full_seq,
                           query_row.sgrna_id])
        if flag:
            df_out = df_row
            flag = False
        else:
            df_out = np.vstack((df_out, df_row))
    df_out = pd.DataFrame(df_out)
    df_out.columns = ['gene_symbol', 'refseq_id', 'exon_id', 'chrom',
                      'start', 'end', 'raw_sequence', 'pam_type',
                      'cutting_site_type', 'cutting_site', 'sgrna_seq',
                      'sgrna_full_seq', 'sgrna_id']
    return df_out


def coordinate_sgrna(df, my_up, my_down, my_sgrna_len):
    for i in range(df.shape[0]):
        record = df.iloc[i, :]
        if record['pam_type'] in ['NGG', 'NAG']:
            old_start = int(record['start'])
            new_start = old_start - (my_sgrna_len - SGRNA_LEN)
            old_raw_seq = record['raw_sequence']
            new_raw_seq = old_raw_seq[-my_sgrna_len:]
            old_sgrna_full_seq = record['sgrna_full_seq']

            start_index = SGRNA_LEN + UPSTREAM - my_sgrna_len - my_up
            end_index = DOWNSTREAM - my_down
            if end_index == 0:
                end_index = None
            new_sgrna_full_seq = old_sgrna_full_seq[start_index:-end_index]

            df.loc[i, 'start'] = new_start
            df.loc[i, 'raw_sequence'] = new_raw_seq
            df.loc[i, 'sgrna_seq'] = new_raw_seq
            df.loc[i, 'sgrna_full_seq'] = new_sgrna_full_seq
        else:
            old_end = int(record['end'])
            new_end = old_end - (SGRNA_LEN - my_sgrna_len)
            old_raw_seq = record['raw_sequence']
            new_raw_seq = old_raw_seq[:my_sgrna_len]
            old_sgrna_seq = record['sgrna_seq']
            new_sgrna_seq = old_sgrna_seq[(SGRNA_LEN - my_sgrna_len):]
            old_sgrna_full_seq = record['sgrna_full_seq']

            start_index = DOWNSTREAM - my_down
            end_index = SGRNA_LEN + UPSTREAM - my_sgrna_len - my_up
            if end_index == 0:
                end_index = None
            new_sgrna_full_seq = old_sgrna_full_seq[start_index:-end_index]

            df.loc[i, 'end'] = new_end
            df.loc[i, 'raw_sequence'] = new_raw_seq
            df.loc[i, 'sgrna_seq'] = new_sgrna_seq
            df.loc[i, 'sgrna_full_seq'] = new_sgrna_full_seq
    return df


def reverse_complement(nt_seq):
    return str(Seq(nt_seq).reverse_complement())
