import sqlalchemy


def gene_symbol_to_refseq(
        genes,
        table_name='igenome_ucsc_hg19_refgene',
        uri='postgresql://yinan:123456@localhost/genome_editing'):
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
