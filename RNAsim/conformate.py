from rnastruct import conformate, mult_conformate
#from RNA import fold_compound, svg_rna_plot
#from seqfold import fold, dot_bracket
import requests
import seqfold as fld
import RNA
import pandas as pd
import numpy as np
import ensembl_rest
import json
from datetime import datetime
import multiprocessing as mp
import traceback

import random, os, shutil

def mainone():
    from RNA import fold_compound, svg_rna_plot
    from seqfold import fold, dot_bracket
    seq = ''.join([random.choice('AUCG') for i in range(random.randint(50, 300))])
    _, my_pred, _, _ = conformate(seq, (3/2)*0.34, nt=0.34, plot=False, processing=False, log=False)
    (vienna_pred, _) = fold_compound(seq).mfe()
    seq_pred = dot_bracket(seq, fold(seq))
    what = mult_conformate(seq)
    print(what)
    print('[INFO] : Generating RNA Plots')
    svg_rna_plot(seq, my_pred, 'my_pred.svg')
    svg_rna_plot(seq, vienna_pred, 'vienna_pred.svg')
    svg_rna_plot(seq, seq_pred, 'seq_pred.svg')

def obey_rls(it, i, j):
    (k, l) = sorted(it)
    (i, j) = sorted((i, j))
    #print(f'i = {i}, j = {j}, k = {k}, l = {l}')
    ret =  (
        (not i in it) and (not j in it) and
        ((k < i and j > l) or
         (i < k and l > j) or
         (j < k and j < l) or
         (i > k and i > l))
    )
    ret2 = (
        (not i in it) and (not j in it) and
        not ((k < i and i < l and l < j) or
         (i < k and k < j and j < l))
    )
    #print(ret, ret2)
    return ret2

def get_plist_from_probs(probs, cutoff=0):
    probs_dict = {}
    for i in range(probs.shape[0]):
        for j in range(probs.shape[1]):
            if probs[i, j] < cutoff:
                continue
            probs_dict[(i, j)] = probs[i, j]
    probs_dict = dict(sorted(probs_dict.items(), key=lambda it: it[1], \
                             reverse=True))
    eps = []
    while len(probs_dict.values()) != 0:
        (i_max, j_max), prob = list(probs_dict.items())[0]
        """
        probs_dict = dict(sorted([it for it in probs_dict.items() \
                                  if not i_max in it[0] and not j_max in it[0]],
                                  key=lambda it: it[1], reverse=True))
        """
        probs_dict = dict(sorted([it for it in probs_dict.items() \
                                  if obey_rls(it[0], i_max, j_max)],
                                  key=lambda it: it[1], reverse=True))
        if prob != 0:
            eps.append(RNA.ep(i_max, j_max, prob, 0))
    """
    for (i, j), prob in probs_dict.items():
        eps.append(RNA.ep(i, j, prob, 0))
    """
    return RNA.ElemProbVector(eps)


def get_avg_conf(seq, **kws):
    matrix = mult_conformate(seq, **kws)
    plist = get_plist_from_probs(matrix)
    logger = kws['logger'] if 'logger' in kws else False
    wr_mat = kws['wr_mat'] if 'wr_mat' in kws else False
    path = os.path.dirname(kws['path']) if 'path' in kws else '.'
    basename = os.path.basename(kws['path'])[:-4] if 'path' in kws else 'avg_conf'
    if wr_mat:
        pd.DataFrame(matrix).to_csv(path + f'/{basename}-matrix.txt', sep='\t')
    if logger:
        for ep in plist:
            print(ep)
    return RNA.db_from_plist(plist, len(seq))

def maintwo(seq, trans, path='.'):
    f = open(path + '/dot_bracket.txt', 'a')
    """
    print('[INFO][{time}]:  Calculating get_avg_conf(trans={transcript}, iter=10)'.format(
        time=datetime.now().strftime("%H:%M:%S"),
        transcript=trans
    ))
    db_iter = get_avg_conf(seq, iter=10, path=path + f'/{trans}-db_iter.png', \
                           wr_mat=True, plot=True)
    f.write(db_iter + '\t' + f'{trans}-' + str('db_iter') + '\t' + \
                    str(len([it for it in db_iter if it != '.'])) + \
                    '\t' + str(len(db_iter)) + '\n')
    print('[INFO][{time}]:  Calculating get_avg_conf(trans={transcript}, iter=1)'.format(
        time=datetime.now().strftime("%H:%M:%S"),
        transcript=trans
    ))
    db_avg = get_avg_conf(seq, iter=1, path=path + f'/{trans}-db_avg.png', \
                          wr_mat=True, plot=True)
    f.write(db_avg + '\t' + f'{trans}-' + str('db_avg') + '\t' + \
                    str(len([it for it in db_avg if it != '.'])) + \
                    '\t' + str(len(db_avg)) + '\n')
    print('[INFO][{time}]:  Calculating conformate(trans={transcript})'.format(
        time=datetime.now().strftime("%H:%M:%S"),
        transcript=trans
    ))
    _, db_one, _, _ = conformate(seq, lp=(3/2)*0.34, plot=True, \
                                 path=path + f'/{trans}-db_one.png')
    f.write(db_one + '\t' + f'{trans}-' + str('db_one') + '\t' + \
                    str(len([it for it in db_one if it != '.'])) + \
                    '\t' + str(len(db_one)) + '\n')
    print('[INFO][{time}]:  Calculating vrna_fold(trans={transcript})'.format(
        time=datetime.now().strftime("%H:%M:%S"),
        transcript=trans
    ))
    db_vrna, _ = RNA.fold_compound(seq).mfe()
    RNA.svg_rna_plot(seq, db_vrna, path + f'/{trans}-db_vrna.svg')
    f.write(db_vrna + '\t' + f'{trans}-' + str('db_vrna') + '\t' + \
                    str(len([it for it in db_vrna if it != '.'])) + \
                    '\t' + str(len(db_vrna)) + '\n')
    """
    print('[INFO][{time}]:  Calculating sqf_fold(trans={transcript})'.format(
        time=datetime.now().strftime("%H:%M:%S"),
        transcript=trans
    ))
    db_sqf = fld.dot_bracket(seq, fld.fold(seq))
    f.write(db_sqf + '\t' + f'{trans}-' + str('db_sqf') + '\t' + \
                    str(len([it for it in db_sqf if it != '.'])) + \
                    '\t' + str(len(db_sqf)) + '\n')
    RNA.svg_rna_plot(seq, db_sqf, path + f'/{trans}-db_sqf.svg')
    """
    mydict = {'db_vrna': db_vrna, 'db_sqf': db_sqf, \
              'db_iter': db_iter, 'db_avg': db_avg, 'db_one': db_one}
    svg = ['db_vrna', 'db_sqf']
    with open(path + '/dot_bracket.txt', 'a') as f:
        for key, i in mydict.items():
            f.write(i + '\t' + f'{trans}-' + str(key) + '\t' + \
                    str(len([it for it in i if it != '.'])) + \
                    '\t' + str(len(i)) + '\n')
            if key in svg:
                RNA.svg_rna_plot(seq, i, path + f'/{trans}-{key}.svg')
    """
    f.close()

def get_sequence(transcript):
    r = requests.get(
        f'https://rest.ensembl.org/sequence/id/{transcript}'\
            '?content-type=text/plain&type=cdna')
    return r.text.replace('T', 'U')

def mp_maintwo(arg_tpl):
    seq, trans, path = arg_tpl
    print('[INFO][{time}]:  Initializing analysis of {transcript}'.format(
                time=datetime.now().strftime("%H:%M:%S"),
                transcript=trans
    ))
    try:
        maintwo(seq, trans, path=path)
        print('[INFO][{time}]:  Completed analysis of {transcript}'.format(
                    time=datetime.now().strftime("%H:%M:%S"),
                    transcript=trans
        ))
    except Exception as e:
        print('[INFO][{time}]:  Error in analysis of {transcript}\n{trace}'.format(
                time=datetime.now().strftime("%H:%M:%S"),
                transcript=trans,
                trace=traceback.format_exc()
        ))
    return 0

def main(path='./results', multi=False):
    if not os.path.exists(path):
        os.mkdir(path)
    with open('mytranscripts.txt', 'r') as f:
        transcripts = f.read().split('\n')
    if multi:
        parallel = []
    for trans in transcripts:
        if not multi:
            print('[INFO][{time}]:  Initializing analysis of {transcript}'.format(
                time=datetime.now().strftime("%H:%M:%S"),
                transcript=trans
            ))
        seq = get_sequence(trans)
        if not os.path.exists(path + f'/{trans}'):
            os.mkdir(path + f'/{trans}')
        with open(path + f'/{trans}/basesequence.txt', 'w') as f:
            f.write(seq)
        json.dump(ensembl_rest.lookup(id=trans), \
                  open(path + f'/{trans}/metadata.json','w'),
                  indent=4)
        if multi:
            parallel.append((seq, trans, path + f'/{trans}'))
        if not multi:
            maintwo(seq, trans, path=path + f'/{trans}')
            print('[INFO][{time}]:  Completed analysis of {transcript}'.format(
                time=datetime.now().strftime("%H:%M:%S"),
                transcript=trans
            ))
    if multi:
        with mp.Pool() as pool:
            res = pool.map(mp_maintwo, parallel)
    print('[INFO][{time}]:  Finished Execution of Program'.format(
        time=datetime.now().strftime("%H:%M:%S")
    ))


if __name__ == '__main__':
    main(multi=False)