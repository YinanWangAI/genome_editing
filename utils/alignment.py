"""Align sequences to reference genome"""
import os
import pandas as pd
import subprocess
import tempfile

BOWTIE_INDEX_PATH = os.getenv('BOWTIE_INDEX_PATH')


def bowtie2_alignment(seq=None, input_file=None, mode='seq', report_all=True):
    bowtie2_index_path = \
        '/Users/yinan/Research/weilab/reference_genome/igenome/UCSC/hg19/Bowtie2Index/genome'
    alignment_out_path = \
        '/Users/yinan/PycharmProjects/genome_editing/tmp/alignment_tmp.sam'
    if mode == 'seq':
        assert seq is not None, 'Please provide sequence'
        if report_all:
            cmd_prefix = 'bowtie2 -a '
        else:
            cmd_prefix = 'bowtie2 '
        cmd = cmd_prefix + '-c -p 2 --end-to-end --score-min L,-1,-1 -x ' + \
              bowtie2_index_path + ' ' + seq + \
              ' -S ' + alignment_out_path
    elif mode == 'file':
        assert input_file is not None, 'Please provide file'
        cmd = 'bowtie2 -a -r -p 2 -x ' + bowtie2_index_path + ' ' + input_file + \
              ' -S ' + alignment_out_path
    else:
        print('Error: Wrong Mode')
        return 1

    run_code = subprocess.call(cmd, shell=True)
    if run_code == 0:
        alignment_out = pd.read_table(alignment_out_path,
                                      comment='@', header=None)
        colnames = ['read_name', 'sum_flags', 'chrom', 'start_1_base',
                    'mapping_quality', 'cigar', 'ref_mates',
                    'start_1_base_mates',
                    'fragment_len', 'read_sequence', 'read_qualities',
                    'alignment_score',
                    'sub_optim_alignment_score', 'num_ambiguous_base',
                    'num_mismatch',
                    'num_gaps', 'num_gap_ext', 'edit_distance', 'mismatch_ref',
                    'filter_reason']
        # alignment_out.columns = colnames
        return alignment_out
    else:
        print('Error in Bowtie2')
        return 1


def bowtie_alignment(seq=None, input_file=None, mode='seq', report_all=True,
                     bowtie_index_path=BOWTIE_INDEX_PATH,
                     off_targets=2):
    temp_file = tempfile.NamedTemporaryFile()
    alignment_out_path = temp_file.name

    if report_all:
        cmd_prefix = 'bowtie -a '
    else:
        cmd_prefix = 'bowtie '

    if mode == 'seq':
        assert seq is not None, 'Please provide sequence'
        cmd = cmd_prefix + '-p 2 -n ' + str(off_targets) + ' -l 23 ' + \
              bowtie_index_path + ' -c ' + seq + \
              ' -S ' + alignment_out_path
    elif mode == 'file':  # TODO: not supported by now
        assert input_file is not None, 'Please provide file'
        cmd = cmd_prefix + '-p 2 -n ' + str(off_targets) + ' -l 23 ' + \
              bowtie_index_path + ' ' + seq + \
              ' -S ' + alignment_out_path
    else:
        print('Error: Wrong Mode')
        return 1

    run_code = subprocess.call(cmd, shell=True)
    if run_code == 0:
        alignment_out = pd.read_table(alignment_out_path,
                                      comment='@', header=None)
        # colnames = ['read_name', 'sum_flags', 'chrom', 'start_1_base',
        #             'mapping_quality', 'cigar', 'ref_mates',
        #             'start_1_base_mates',
        #             'fragment_len', 'read_sequence', 'read_qualities',
        #             'alignment_score',
        #             'sub_optim_alignment_score', 'num_ambiguous_base',
        #             'num_mismatch',
        #             'num_gaps', 'num_gap_ext', 'edit_distance', 'mismatch_ref',
        #             'filter_reason']
        # # alignment_out.columns = colnames
        return alignment_out
    else:
        print('Error in Bowtie')
        return 1
