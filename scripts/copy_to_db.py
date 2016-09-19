"""Copy data from files to database
sqlalchemy URI: postgresql://wangyn@localhost/genome_editing (knight)
postgresql://yinan:123456@localhost/genome_editing (local)
"""
import sqlalchemy
import genome_editing.design_sgRNA

# copy whole genome sgRNAs into database, only store the upstream 20bp
whole_genome_sgrna_mat = whole_genome_sgrna()
with open('./whole_genome_sgrna_mat.pkl', 'rb') as f:
    whole_genome_sgrna_mat = pickle.load(f)
engine = sqlalchemy.create_engine(
    'postgresql://wangyn@localhost/genome_editing')
whole_genome_sgrna_gg = pd.DataFrame(whole_genome_sgrna_mat[0])
whole_genome_sgrna_gg.to_sql('grch38_wg_sgrna_ngg', engine,
                             index=None, if_exists='replace', chunksize=10000)
