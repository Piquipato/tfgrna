from textwrap import wrap as splitk
import ensembl_rest
import pandas as pd
import os, shutil
import requests

def get_transcript(transcript):
    r = requests.get(
        f'https://rest.ensembl.org/sequence/id/{transcript}'\
            '?content-type=text/plain&type=cdna')
    return r.text.replace('T', 'U')

def load_dot_bracket(file):
    return pd.read_csv(file, sep='\t', \
            header=None, index_col=None)

beauty = lambda s, x: '\n'.join(splitk(s, x))

def main(path="./dot_bracket"):
    for file in [it for it in os.listdir(path) if '.txt' in it]:
        name = file[:-4]
        db = load_dot_bracket(path + '/' + file)
        transcript = list(db[1].values)[0].split('-')[0]
        sequence = get_transcript(transcript)
        wrt_str = '[SEQUENCE]\n'
        wrt_str += beauty(sequence, 80)
        wrt_str += '\n'
        dots = list(db[0].values)
        where = [item.split('-')[1] for item in list(db[1].values)]
        for i in range(len(dots)):
            wrt_str += '[{where}]\n'.format(where=where[i].upper())
            wrt_str += beauty(dots[i], 80)
            wrt_str += '\n'
        with open(path + f'/dot_bracket-{name}.txt', 'w') as f:
            f.write(wrt_str)


if __name__ == '__main__':
    main()