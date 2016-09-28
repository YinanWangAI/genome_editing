"""Design sgRNAs for CRISPR/Cas9 gene editing"""
from Bio.Seq import Seq
import numpy as np
import pandas as pd
import regex
import sqlalchemy
from .rs2 import compute_rs2


class Designer:
    """Design sgRNAs for a target"""

    def __init__(self, entrez_id, sgrna_upstream=4, sgrna_downstream=3,
                 sgrna_length=20, flank=30, overlapped=True, filter_tttt=True):
        """Init

        Args:
            entrez_id: the Entrez ID of target gene
            sgrna_upstream: the length of sgRNA upstream
            sgrna_downstream: the length of sgRNA downstream
            flank: the length of flank sequence
            overlapped: whether the sgRNAs could be overlapped
        """
        self.target_gene = Gene(entrez_id)
        self.target_gene.get_sequence(flank)
        self.sgrna_upstream = sgrna_upstream
        self.sgrna_downstream = sgrna_downstream
        self.sgrna_length = sgrna_length
        self.flank = flank
        self.sgrnas = []
        self.overlapped = overlapped
        self.filter_tttt = filter_tttt

    def get_sgrnas(self, pams=['NGG', 'NAG']):
        """Get sgRNAs with PAM NGG and NAG

        Args:
            pams: the PAM to design

        Returns:
            None. The results are stored in self.sgrnas
        """
        exon_num = self.target_gene.exons.shape[0]
        for i in range(exon_num):
            exon_seq = self.target_gene.exons.seq[i]
            exon_start = self.target_gene.exons.start[i]
            sgrnas = []
            for pam in pams:
                pam_pattern = self._get_sgrna_pattern(pam)
                pam_pattern_rc = self._get_sgrna_pattern(
                    pam, reverse_complement=True)
                sgrnas += self._design_sgrna(exon_seq, pam_pattern, pam)
                sgrnas += self._design_sgrna(exon_seq, pam_pattern_rc,
                                             self._reverse_complement(pam),
                                             reverse_complement=True)
            for sgrna in sgrnas:
                sgrna.chrom = self.target_gene.chrom
                sgrna.entrez_id = self.target_gene.entrez_id
                sgrna.exon_id = i
                sgrna.start += exon_start - self.flank
                sgrna.end += exon_start - self.flank
                if sgrna.rc:
                    sgrna.cutting_site = sgrna.start + 2.5
                else:
                    sgrna.cutting_site = sgrna.end - 2.5
                self.sgrnas.append(sgrna)

    def _get_sgrna_pattern(self, pam, reverse_complement=False):
        if reverse_complement:
            pam = self._reverse_complement(pam)
            pam_pattern = ''
            for char in pam:
                if char == 'N':
                    pam_pattern += '\w'
                else:
                    pam_pattern += char
            sgrna_pattern = '\w{%d}%s\w{%d}\w{%d}' % (
                self.sgrna_downstream, pam_pattern,
                self.sgrna_length, self.sgrna_upstream
            )
        else:
            pam_pattern = ''
            for char in pam:
                if char == 'N':
                    pam_pattern += '\w'
                else:
                    pam_pattern += char
            sgrna_pattern = '\w{%d}\w{%d}%s\w{%d}' % (
                self.sgrna_upstream, self.sgrna_length,
                pam_pattern, self.sgrna_downstream
            )
        return sgrna_pattern

    def _design_sgrna(self, seq, pam_pattern, pam, reverse_complement=False):
        sgrna_match = regex.finditer(pam_pattern, seq,
                                     overlapped=self.overlapped)
        sgrnas = []
        for sgrna in sgrna_match:
            full_seq = sgrna.group()
            if reverse_complement:
                sgrna_seq = full_seq[(self.sgrna_downstream +
                                      len(pam)):-self.sgrna_upstream]
            else:
                sgrna_seq = full_seq[self.sgrna_upstream:(self.sgrna_upstream +
                                                          self.sgrna_length)]
            if self.filter_tttt:
                if reverse_complement:
                    filter_pattern = 'AAAA'
                else:
                    filter_pattern = 'TTTT'
                if sgrna_seq.find(filter_pattern) != -1:
                    continue
            if reverse_complement:
                sgrna_start = sgrna.start() + len(pam) + self.sgrna_downstream
                sgrna_end = sgrna_start + self.sgrna_length - 1
                sgrna_cutting_site = sgrna.start() + 2.5
            else:
                sgrna_start = sgrna.start() + self.sgrna_upstream
                sgrna_end = sgrna_start + self.sgrna_length - 1
                sgrna_cutting_site = sgrna_end - 2.5
            if (sgrna_cutting_site < self.flank) or \
                    (sgrna_cutting_site >= (len(seq) - self.flank)):
                sgrna_type = 'splicing site'
            else:
                sgrna_type = 'coding region'
            sgrnas.append(SgRNA(sequence=sgrna_seq,
                                cutting_site_type=sgrna_type,
                                start=sgrna_start, end=sgrna_end,
                                pam_type=pam,
                                full_seq=full_seq,
                                rc=reverse_complement))
        return sgrnas

    def _reverse_complement(self, sgrna_seq):
        """Get reverse complement of sequence

        Args:
            sgrna_seq: sequence

        Returns:
            str, the reverse complement of input sequence
        """
        return str(Seq(sgrna_seq).reverse_complement())

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
            if sgrna.rc:
                seq = sgrna.reverse_complement()
            else:
                seq = sgrna.sequence
            df_row = [sgrna.entrez_id, sgrna.exon_id, sgrna.chrom, sgrna.start,
                      sgrna.end, sgrna.sequence, sgrna.pam_type,
                      sgrna.cutting_site_type, sgrna.cutting_site, seq,
                      sgrna.full_seq]
            if flag:
                df = np.asarray(df_row)
                flag = False
            else:
                df = np.vstack((df, df_row))
        df = pd.DataFrame(df)
        df.columns = ['entrez_id', 'exon_id', 'chrom', 'start', 'end',
                      'raw_sequence', 'pam_type', 'cutting_site_type',
                      'cutting_site', 'sgrna_seq', 'sgrna_full_seq']
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
                 exon_id=None, cutting_site=None, full_seq=None,
                 aa_cut=None, per_peptide=None, rs2_score=None,
                 rc=None):
        """Init

        Args:
            sequence: sgRNA sequence
            type: the type of sgRNA,
            (3'UTR, 5'UTR, splicing site, coding region)
            gene_id: the gene to be targeted
            chrom: location
            start: start position, 1 based
            end: end position, 1 based
            rc: if in reverse_complement strand
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
        self.full_seq = full_seq
        self.aa_cut = aa_cut
        self.per_peptide = per_peptide
        self.rc = rc
        # compute rs2 score
        # self.rs2_score = rs2_score
        # if rs2_score is not None:
        #     self.rs2_score = rs2_score
        # else:
        #     self.rs2_score = compute_rs2(self.sequence, self.aa_cut,
        #                                  self.per_peptide)

    def get_gc_content(self):
        if self.sequence is None:
            return None
        else:
            gc_content = (self.sequence.count('G') +
                          self.sequence.count('C')) / len(self.sequence)
            return gc_content

    def reverse_complement(self):
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

    def get_rs2_score(self):
        assert self.pam_type == 'NGG', 'Only support NGG'
        assert len(self.full_seq) == 30, 'Have to provide 30mers'
        return compute_rs2(self.full_seq, self.aa_cut, self.per_peptide)
