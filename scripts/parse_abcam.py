"""Extract positive control and Entrez ID from abcam web page"""
from bs4 import BeautifulSoup
import regex


def get_positive_control(html_path):
    """Extract positive controls

    Args:
        html_path: the path of html file

    Returns:
        positive controls
    """
    soup = BeautifulSoup(open(html_path), 'lxml')
    li_items = soup.find_all('li')
    flag = False
    for li_item in li_items:
        text = li_item.get_text()
        if flag:
            return text
        elif 'Positive control' in text:
            flag = True
    return None


def get_entrez_id(html_path):
    """Extract Entrez ID

    Args:
        html_path: the path of html file

    Returns:
        Entrez ID
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
    html_path = '/Users/yinan/Desktop/53bp1-antibody-ab36823.html'
    print(get_positive_control(html_path))
    print(get_entrez_id(html_path))
