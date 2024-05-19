import pandas as pd
import numpy as np
import os
import json
import seqfold as fld
import langevin as lg
import RNA


def get_trans_data(path):
    trs = {}
    with open(path + '/basesequence.txt', 'r') as f:
        trs['seq'] = f.read()
    trs['dot_bracket'] = {}
    with open(path + '/dot_bracket-db_rndc.txt', 'r') as g:
        trs['dot_bracket']['db_rndc'] = ''.join(g.read().split('[db_rndc]')[-1].split('\n'))
    with open(path + '/dot_bracket.txt', 'r') as b:
        for line in [it for it in b.read().split('\n') if it != '']:
            trs['dot_bracket'][line.split('\t')[1].split('-')[-1]] = line.split('\t')[0]
    trs['tbl'] = {}
    trs['tbl']['db_one'] = pd.read_excel(path + f'/{os.path.basename(path)}-db_one-nucleotides.xlsx')
    trs['tbl']['db_iter'] = pd.read_excel(path + f'/{os.path.basename(path)}-db_iter-nucleotides.xlsx')
    trs['name'] = json.load(open(path + '/metadata.json', 'r'))['display_name']
    trs['id'] = os.path.basename(path)
    return trs

def get_results_data(path='./results'):
    results = []
    for dir in os.listdir(path):
        new_path = os.path.join(os.path.abspath(path), dir)
        if os.path.isdir(new_path):
            results.append(get_trans_data(new_path))
    return results

mod = lambda x: np.sqrt(np.dot(x, x))

def calc_g_free(data, nt=0.34):
    R = 1.9872e-3 # kcal/K mol
    T = 37
    rndw = ['db_one', 'db_iter']
    dynp = ['db_rndc', 'db_sqf', 'db_vrna']
    N = len(data['seq'])
    dg = {}
    df = []
    for db in data['dot_bracket'].keys():
        binds = len([it for it in data['dot_bracket'][db] if it != '.'])
        if db in dynp:
            g = RNA.eval_structure_simple(data['seq'], data['dot_bracket'][db])
            dg[db] = g / (R*T)
            df.append([db, dg[db], binds, len(data['seq'])])
        if db in rndw:
            row = data['tbl'][db].iloc[-1, :]
            row = [float(it) for it in row['Position'][1:-1].split(' ') \
                   if it != '']
            pos = np.array(row, dtype=np.double)
            dg[db] = - N * lg.lnZ(lg.inv_lange(mod(pos) / (N*nt)))
            df.append([db, dg[db], binds, len(data['seq'])])
            g = RNA.eval_structure_simple(data['seq'], data['dot_bracket'][db])
            dg[db+'_rna'] = g / (R*T)
            df.append([db+'_rna', dg[db+'_rna'], binds, len(data['seq'])])
    data['dg'] = dg
    data['df'] = pd.DataFrame(df, \
                              columns=pd.MultiIndex.from_tuples([(data['name'], x) for x in ['model', 'dg', 'binds', 'length']]))
    return data



if __name__ == '__main__':
    res = get_results_data()
    dfs = pd.concat([calc_g_free(r)['df'] for r in res], axis=1)
    with pd.ExcelWriter('freeg.xlsx', mode='a', if_sheet_exists='new') as writer:
        dfs.to_excel(writer, sheet_name='Free Energy')