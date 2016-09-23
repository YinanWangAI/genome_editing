"""Extract positive control and Entrez ID from abcam web page"""
import argparse
from bs4 import BeautifulSoup
import numpy as np
import os
import pandas as pd
import regex


def batch_parse(input_dir, output_path='./output.csv'):
    """Parse abcam html files

    Args:
        input_dir: the dir storing files
        output_path: the path of output

    Returns:
        None
    """
    file_path = os.listdir(input_dir)
    flag = True
    for i in range(len(file_path)):
        if i % 5000 == 0:
            print('----------{}%'.format(round(i / len(file_path) * 100, 2)))
        single_path = os.path.join(input_dir, file_path[i])
        positive_controls = get_positive_control(single_path)
        gene_id = get_entrez_id(single_path)
        for cell_line in positive_controls:
            if flag:
                flag = False
                output = np.array([gene_id, cell_line])
            else:
                output = np.vstack((output, np.array([gene_id, cell_line])))
    output = pd.DataFrame(output)
    output.columns = ['gene_id', 'positive_control']
    output.to_csv(output_path, index=None)


def get_positive_control(html_path):
    """Extract positive controls

    Args:
        html_path: file path

    Returns:
        positive_control, a list of positive controls
    """
    soup = BeautifulSoup(open(html_path), 'lxml')
    li_items = soup.find_all('li')
    flag = False
    for li_item in li_items:
        text = li_item.get_text()
        if flag:
            positive_control = []
            # split by ','
            cell_lines = text.split(',')
            for i in range(len(cell_lines)):
                if cell_lines[i][0] == ' ':
                    cell_lines[i] = cell_lines[i][1:]
            # split by 'and'
            for cell_line in cell_lines:
                cell_line = cell_line.split('and')
                for sub_cell_line in cell_line:
                    if sub_cell_line[0] == ' ':
                        sub_cell_line = sub_cell_line[1:]
                    if (sub_cell_line[-1] == '.') | (sub_cell_line[-1] == ' '):
                        sub_cell_line = sub_cell_line[:-1]
                    positive_control.append(sub_cell_line)
            return positive_control
        elif 'Positive control' in text:
            flag = True
    return None


def get_entrez_id(html_path):
    """Extract Entrez ID

    Args:
        html_path: file path

    Returns:
        gene id
    """
    soup = BeautifulSoup(open(html_path), 'lxml')
    li_items = soup.find_all('li')
    pattern = '^EntrezGene:(\d*)Human$'
    for li_item in li_items:
        text = li_item.get_text()
        text = text.replace('\n', '')
        text = text.replace(' ', '')
        m = regex.match(pattern, text)
        if m is not None:
            return m.group(1)
    return None


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='batch_parse parameters')
    parser.add_argument('path', metavar='N', type=str, nargs='+',
                        help='batch_parse parameters: input and output path')
    args = parser.parse_args()
    # print(args.path)
    batch_parse(args.path[0], args.path[1])
