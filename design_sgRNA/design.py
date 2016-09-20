"""Design sgRNAs for CRISPR/Cas9 gene editing"""
from Bio.Seq import Seq
import numpy as np
import pandas as pd
import regex
import sqlalchemy


class Designer:
    """Design sgRNAs for a target"""

    def __init__(self, entrez_id, sgrna_upstream=20, sgrna_downstream=3,
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
        """Get sgRNAs with PAM NGG and NAG

        Args:
            filter_tttt: whether filter sgRNAs containing TTTT

        Returns:
            None. The results are stored in self.sgrnas
        """
        exon_num = self.target_gene.exons.shape[0]
        for i in range(exon_num):
            exon_seq = self.target_gene.exons.seq[i]
            exon_start = self.target_gene.exons.start[i]
            sgrnas_ngg = self._design_sgrna_ngg(exon_seq, filter_tttt)
            sgrnas_nag = self._design_sgrna_nag(exon_seq, filter_tttt)
            for sgrna in (sgrnas_ngg + sgrnas_nag):
                sgrna.chrom = self.target_gene.chrom
                sgrna.entrez_id = self.target_gene.entrez_id
                sgrna.exon_id = i
                sgrna.start += exon_start - self.flank
                sgrna.end += exon_start - self.flank
                if sgrna.pam_type == 'NGG' or sgrna.pam_type == 'NAG':
                    sgrna.cutting_site = sgrna.end - 2.5
                elif sgrna.pam_type == 'CCN' or sgrna.pam_type == 'CTN':
                    sgrna.cutting_site = sgrna.start + 2.5
                self.sgrnas.append(sgrna)

    def _reverse_complement(self, sgrna_seq):
        """Get reverse complement of sequence

        Args:
            sgrna_seq: sequence

        Returns:
            str, the reverse complement of input sequence
        """
        return str(Seq(sgrna_seq).reverse_complement())

    def _design_sgrna_ngg(self, seq, filter_tttt=True):
        return self._ngg(seq, filter_tttt) + self._ccn(seq, filter_tttt)

    def _design_sgrna_nag(self, seq, filter_tttt=True):
        return self._nag(seq, filter_tttt) + self._ctn(seq, filter_tttt)

    def _ngg(self, seq, filter_tttt=True):
        """Design sgRNAs with sequence *NGG*

        Args:
            seq: the sequence to be targeted
            filter_tttt: whether filter sgRNAs containing TTTT

        Returns:
            sgrnas, a list of SgRNA object
        """

        sgrna_match = regex.finditer(
            '\w{%d}\wGG\w{%d}' % (self.sgrna_upstream, self.sgrna_downstream),
            seq, overlapped=True)
        sgrnas = []
        for sgrna in sgrna_match:
            sgrna_seq = sgrna.group()[:self.sgrna_upstream]
            if filter_tttt:
                if sgrna_seq.find('TTTT') != -1:
                    continue
            sgrna_start = sgrna.start()
            sgrna_end = sgrna_start + self.sgrna_upstream - 1
            sgrna_cutting_site = sgrna_end - 2.5
            if (sgrna_cutting_site < self.flank) or \
                    (sgrna_end >= (len(seq) - self.flank)):
                sgrna_type = 'splicing site'
            else:
                sgrna_type = 'coding region'
            sgrnas.append(SgRNA(sequence=sgrna_seq,
                                cutting_site_type=sgrna_type,
                                start=sgrna_start, end=sgrna_end,
                                pam_type='NGG'))
        return sgrnas

    def _ccn(self, seq, filter_tttt=True):
        """Design sgRNAs with sequence *CCN*

                Args:
                    seq: the sequence to be targeted
                    filter_tttt: whether filter sgRNAs containing TTTT

                Returns:
                    sgrnas, a list of SgRNA object
                """
        sgrna_match = regex.finditer(
            '\w{%d}CC\w\w{%d}' % (self.sgrna_downstream, self.sgrna_upstream),
            seq, overlapped=True)
        sgrnas = []
        for sgrna in sgrna_match:
            sgrna_seq = sgrna.group()[-self.sgrna_upstream:]
            if filter_tttt:
                if sgrna_seq.find('AAAA') != -1:
                    continue
            sgrna_start = sgrna.end() - self.sgrna_upstream
            sgrna_end = sgrna.end() - 1
            sgrna_cutting_site = sgrna_start + 2.5
            if (sgrna_cutting_site < self.flank) or \
                    (sgrna_cutting_site >= (len(seq) - self.flank)):
                sgrna_type = 'splicing site'
            else:
                sgrna_type = 'coding region'
            sgrnas.append(SgRNA(sequence=sgrna_seq,
                                cutting_site_type=sgrna_type,
                                start=sgrna_start, end=sgrna_end,
                                pam_type='CCN'))
        return sgrnas

    def _nag(self, seq, filter_tttt=True):
        """Design sgRNAs with sequence *NAG*

        Args:
            seq: the sequence to be targeted
            filter_tttt: whether filter sgRNAs containing TTTT

        Returns:
            sgrnas, a list of SgRNA object
        """

        sgrna_match = regex.finditer(
            '\w{%d}\wAG\w{%d}' % (self.sgrna_upstream, self.sgrna_downstream),
            seq, overlapped=True)
        sgrnas = []
        for sgrna in sgrna_match:
            sgrna_seq = sgrna.group()[:self.sgrna_upstream]
            if filter_tttt:
                if sgrna_seq.find('TTTT') != -1:
                    continue
            sgrna_start = sgrna.start()
            sgrna_end = sgrna_start + self.sgrna_upstream - 1
            sgrna_cutting_site = sgrna_end - 2.5
            if (sgrna_cutting_site < self.flank) or \
                    (sgrna_end >= (len(seq) - self.flank)):
                sgrna_type = 'splicing site'
            else:
                sgrna_type = 'coding region'
            sgrnas.append(SgRNA(sequence=sgrna_seq,
                                cutting_site_type=sgrna_type,
                                start=sgrna_start, end=sgrna_end,
                                pam_type='NAG'))
        return sgrnas

    def _ctn(self, seq, filter_tttt=True):
        """Design sgRNAs with sequence *CTN*

                Args:
                    seq: the sequence to be targeted
                    filter_tttt: whether filter sgRNAs containing TTTT

                Returns:
                    sgrnas, a list of SgRNA object
                """
        sgrna_match = regex.finditer(
            '\w{%d}CT\w\w{%d}' % (self.sgrna_downstream, self.sgrna_upstream),
            seq, overlapped=True)
        sgrnas = []
        for sgrna in sgrna_match:
            sgrna_seq = sgrna.group()[-self.sgrna_upstream:]
            if filter_tttt:
                if sgrna_seq.find('AAAA') != -1:
                    continue
            sgrna_start = sgrna.end() - self.sgrna_upstream
            sgrna_end = sgrna.end() - 1
            sgrna_cutting_site = sgrna_start + 2.5
            if (sgrna_cutting_site < self.flank) or \
                    (sgrna_cutting_site >= (len(seq) - self.flank)):
                sgrna_type = 'splicing site'
            else:
                sgrna_type = 'coding region'
            sgrnas.append(SgRNA(sequence=sgrna_seq,
                                cutting_site_type=sgrna_type,
                                start=sgrna_start, end=sgrna_end,
                                pam_type='CTN'))
        return sgrnas

    def print(self):
        """Print the sgRNAs

        Returns:
            None
        """
        for sgrna in self.sgrnas:
            sgrna.print()

    def output(self):
        """Output sgRNAs in a pandas DataFrame

        Returns:
            pd.DataFrame
        """
        flag = True
        for sgrna in self.sgrnas:
            df_row = [sgrna.entrez_id, sgrna.exon_id, sgrna.chrom, sgrna.start,
                      sgrna.end, sgrna.sequence, sgrna.pam_type,
                      sgrna.cutting_site_type, sgrna.cutting_site]
            if flag:
                df = np.asarray(df_row)
                flag = False
            else:
                df = np.vstack((df, df_row))
        df = pd.DataFrame(df)
        df.columns = ['entrez_id', 'exon_id', 'chrom', 'start', 'end',
                      'sequence', 'pam_type', 'cutting_site_type',
                      'cutting_site']
        df.loc[:, 'cutting_site'] = df.cutting_site.astype(np.float)
        return df

    def print_cutting_site(self):
        sgrnas_df = self.output()
        cutting_site_coding = sgrnas_df[
            sgrnas_df.cutting_site_type == 'coding region']
        for i in range(self.target_gene.exons.shape[0]):
            exon_start = self.target_gene.exons.start[i]
            exon_seq = self.target_gene.exons.seq[i][self.flank:-self.flank]
            exon_cutting_site = np.floor(
                cutting_site_coding[
                    cutting_site_coding.exon_id == str(i)].cutting_site.values)
            exon_cutting_site = exon_cutting_site - exon_start
            exon_cutting_site.sort()
            exon_with_cutting = ''
            break_start = 0
            for cut_index in exon_cutting_site:
                break_end = int(cut_index) + 1
                exon_with_cutting = exon_with_cutting + \
                                    exon_seq[break_start:break_end] + '    '
                break_start = break_end
            exon_with_cutting = exon_with_cutting + exon_seq[break_start:]
            print('exon id: {}'.format(i))
            print(exon_with_cutting)
            print('\n\n\n\n')


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
        self.exons.index = range(self.exons.shape[0])
        self.chrom = list(set(self.exons.loc[:, 'chrom'].values))
        assert len(self.chrom) == 1, print('Multiple chromosomes')
        self.chrom = self.chrom[0]
        self.gene_symbol = list(set(self.exons.gene_symbol.values))
        self.exons_start = np.min(self.exons.start)
        self.exons_end = np.max(self.exons.end)

    def get_sequence(self, flank):
        """Get exons' sequences with flank

        Args:
            flank: length of flank

        Returns:
            None, the results are stored in self.exons
        """
        self.exons.loc[:, 'seq'] = ''
        table_name = 'grch38_' + self.chrom
        chrom_seq = pd.read_sql(table_name, self.engine).iloc[0, 0]
        for i in range(self.exons.shape[0]):
            start = self.exons.loc[:, 'start'].values[i] - flank
            end = self.exons.loc[:, 'end'].values[i] + flank
            self.exons.loc[i, 'seq'] = chrom_seq[(start - 1):end].upper()


class SgRNA:
    """SgRNA targeting a region of the genome"""

    def __init__(self, sequence=None, pam_type=None, cutting_site_type=None,
                 entrez_id=None, chrom=None, start=None, end=None,
                 exon_id=None, cutting_site=None):
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
        self.pam_type = pam_type
        self.entrez_id = entrez_id
        self.chrom = chrom
        self.start = start
        self.end = end
        self.cutting_site = cutting_site
        self.cutting_site_type = cutting_site_type
        self.exon_id = exon_id

    def get_gc_content(self):
        if self.sequence is None:
            return None
        else:
            gc_content = (self.sequence.count('G') +
                          self.sequence.count('C')) / len(self.sequence)
            return gc_content

    def _reverse_complement(self):
        """Get reverse complement of sequence

        Args:
            sgrna_seq: sequence

        Returns:
            str, the reverse complement of input sequence
        """
        return str(Seq(self.sequence).reverse_complement())

    def print(self):
        """Print the SgRNA object

        Returns:
            None
        """
        print('entrez_id: {}, exon_id: {}, chrom: {}, start: {}, end: {}, '
              'sequence: {}, PAM_type: {}, cutting_site_type: {}, '
              'cutting_site: {}'.format(
            self.entrez_id, self.exon_id, self.chrom, self.start, self.end,
            self.sequence, self.pam_type, self.cutting_site_type,
            self.cutting_site
        ))
