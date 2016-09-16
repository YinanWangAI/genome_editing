"""Design sgRNAs for CRISPR/Cas9 gene editing"""
from Bio.Seq import Seq
import numpy as np
import pandas as pd
import regex
import sqlalchemy


class Designer:
    """Design sgRNAs for a target"""

    def __init__(self, entrez_id, sgrna_upstream=24, sgrna_downstream=3,
                 flank=30):
        """Init

        Args:
            entrez_id: the Entrez ID of target gene
            sgrna_upstream: the length of sgRNA upstream
            sgrna_downstream: the length of sgRNA downstream
            flank: the length of flank sequence
        """
        self.target_gene = Gene(entrez_id)
        self.target_gene.get_sequence(flank)
        self.sgrna_upstream = sgrna_upstream
        self.sgrna_downstream = sgrna_downstream
        self.flank = flank
        self.sgrnas = []

    def get_sgrnas(self, filter_tttt=True):
        """Get sgRNAs

        Args:
            filter_tttt: whether filter sgRNAs containing TTTT

        Returns:
            None. The results are stored in self.sgrnas
        """
        exon_num = self.target_gene.exons.shape[0]
        for i in range(exon_num):
            exon_seq = self.target_gene.exons.seq[i]
            exon_start = self.target_gene.exons.start[i]
            sgrnas_gg = self._design_sgrna_gg(exon_seq, filter_tttt)
            exon_seq_rc = self._reverse_complement(exon_seq)
            sgrnas_cc = self._design_sgrna_cc(exon_seq_rc, filter_tttt)
            for sgrna in (sgrnas_gg + sgrnas_cc):
                sgrna.chrom = self.target_gene.chrom
                sgrna.entrez_id = self.target_gene.entrez_id
                sgrna.exon_id = i
                sgrna.start += exon_start - self.flank
                sgrna.end += exon_start - self.flank
                self.sgrnas.append(sgrna)

    def _reverse_complement(self, sgrna_seq):
        """Get reverse complement of sequence

        Args:
            sgrna_seq: sequence

        Returns:
            str, the reverse complement of input sequence
        """
        return str(Seq(sgrna_seq).reverse_complement())

    def _design_sgrna_gg(self, seq, filter_tttt=True):
        """Design sgRNAs with sequence *NGG*

        Args:
            seq: the sequence to be targeted
            filter_tttt: whether filter sgRNAs containing TTTT

        Returns:
            sgrnas, a list of SgRNA object
        """

        sgrna_match = regex.finditer(
            '\w{%d}\wGG\w{%d}'%(self.sgrna_upstream, self.sgrna_downstream),
            seq)
        sgrnas = []
        for sgrna in sgrna_match:
            sgrna_seq = sgrna.group()
            if filter_tttt:
                if sgrna_seq.find('TTTT') != -1:
                    continue
            sgrna_start = sgrna.start()
            sgrna_end = sgrna.end() - 1
            if (sgrna_start < self.flank) or \
                    (sgrna_end >= (len(seq) - self.flank)):
                sgrna_type = 'splicing site'
            else:
                sgrna_type = 'coding region'
            sgrnas.append(SgRNA(sequence=sgrna_seq, sgrna_type=sgrna_type,
                                start=sgrna_start, end=sgrna_end))
        return sgrnas

    def _design_sgrna_cc(self, seq, filter_tttt=True):
        """Design sgRNAs with sequence *CCN*

                Args:
                    seq: the sequence to be targeted
                    filter_tttt: whether filter sgRNAs containing TTTT

                Returns:
                    sgrnas, a list of SgRNA object
                """
        sgrna_match = regex.finditer(
            '\w{%d}CC\w\w{%d}'%(self.sgrna_downstream, self.sgrna_upstream),
            seq)
        sgrnas = []
        for sgrna in sgrna_match:
            sgrna_seq = sgrna.group()
            if filter_tttt:
                if sgrna_seq.find('TTTT') != -1:
                    continue
            sgrna_start = sgrna.start()
            sgrna_end = sgrna.end() - 1
            if (sgrna_start < self.flank) or \
                    (sgrna_end >= (len(seq) - self.flank)):
                sgrna_type = 'splicing site'
            else:
                sgrna_type = 'coding region'
            sgrnas.append(SgRNA(sequence=sgrna_seq, sgrna_type=sgrna_type,
                                start=sgrna_start, end=sgrna_end))
        return sgrnas

    def print(self):
        """Print the sgRNAs

        Returns:
            None
        """
        for sgrna in self.sgrnas:
            sgrna.print()


class Gene:
    """The gene to be edited"""

    def __init__(self, entrez_id,
                 table_name='grch38_exon',
                 engine=sqlalchemy.create_engine(
                     'postgresql://yinan:123456@localhost/genome_editing')
                 ):
        """Init

        Args:
            entrez_id: the ENTREZ ID of the gene
            table_name: the table in the database storing exon information
            engine: sqlalchemy engine
        """
        self.entrez_id = entrez_id
        query = "SELECT * FROM {} WHERE entrez_id='{}'".format(table_name,
                                                               str(entrez_id))
        self.engine = engine
        self.exons = pd.read_sql_query(
            query, self.engine).drop_duplicates().sort_values('start')
        self.chrom = list(set(self.exons.loc[:, 'chrom'].values))
        assert len(self.chrom) == 1, print('Multiple chromosomes')
        self.chrom = self.chrom[0]
        self.gene_symbol = list(set(self.exons.gene_symbol.values))
        self.exons_start = np.min(self.exons.start)
        self.exons_end = np.max(self.exons.end)

    def get_sequence(self, flank=30):
        """Get exons' sequences with flank

        Args:
            flank: length of flank

        Returns:
            None, the results are stored in self.exons
        """
        self.exons.loc[:, 'seq'] = ''
        self.flank = flank
        table_name = 'grch38_' + self.chrom
        chrom_seq = pd.read_sql(table_name, self.engine).iloc[0, 0]
        for i in range(self.exons.shape[0]):
            start = self.exons.loc[:, 'start'][i] - flank
            end = self.exons.loc[:, 'end'][i] + flank
            self.exons.loc[i:, 'seq'] = chrom_seq[(start - 1):end].upper()


class SgRNA:
    """SgRNA targeting a region of the genome"""

    def __init__(self, sequence=None, sgrna_type=None, entrez_id=None,
                 chrom=None, start=None, end=None, exon_id=None):
        """Init

        Args:
            sequence: sgRNA sequence
            type: the type of sgRNA,
            (3'UTR, 5'UTR, splicing site, coding region)
            gene_id: the gene to be targeted
            chrom: location
            start: start position, 1 based
            end: end position, 1 based
        """
        self.sequence = sequence
        self.sgrna_type = sgrna_type
        self.entrez_id = entrez_id
        self.chrom = chrom
        self.start = start
        self.end = end
        self.exon_id = exon_id

    def print(self):
        """Print the SgRNA object

        Returns:
            None
        """
        print('sequence: {}, chrom: {}, type: {}, entrez_id: {}, '
              'start: {}, end: {}, exon_id: {}'.format(
            self.sequence, self.chrom, self.sgrna_type,
            self.entrez_id, self.start, self.end, self.exon_id))
