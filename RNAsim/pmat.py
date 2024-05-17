import numpy as np
import random
import requests
from textwrap import wrap as splitk
import RNA

def p_mat_conf(sequence):
  compl = {
      'A': ['T', 'U'],
      'T': ['A'],
      'U': ['A', 'G'],
      'G': ['C', 'U'],
      'C': ['G']
  }
  n = len(sequence)
  p = np.zeros((n, n))
  idxs = []
  for i in range(n):
    for j in range(n):
      idxs.append((i, j))
      if i != j and sequence[i] in compl[sequence[j]] and abs(i-j) > 3:
        p[i, j] += 1
        p[j, i] += 1
  p /= np.sum(p)
  p_ij = dict(zip(idxs, p.flatten()))
  chosen = []
  while len(list(p_ij.values())) != 0:
    x = np.random.choice(range(len(list(p_ij.keys()))), p=list(p_ij.values()))
    (i, j) = list(p_ij.keys())[x]
    p_ij = {(k, l): val for (y, ((k, l), val)) in enumerate(p_ij.items()) \
            if (
                k != i and k != j and l != i and l != j and \
                not ((k in range(*sorted((i, j))) and not l in range(*sorted((i, j)))) or
                     (not k in range(*sorted((i, j))) and l in range(*sorted((i, j)))))
            )}
    keys = list(p_ij.keys())
    vals = list(p_ij.values())
    if np.sum(vals) == 0:
      break
    vals /= np.sum(vals)
    p_ij = dict(zip(keys, vals))
    chosen.append((i, j))
  dot_bracket = ['.' for k in range(n)]
  for ch in chosen:
    (i, j) = (np.min(ch), np.max(ch))
    dot_bracket[i] = '('
    dot_bracket[j] = ')'
  return ''.join(dot_bracket), chosen

def get_sequence(transcript):
    r = requests.get(
        f'https://rest.ensembl.org/sequence/id/{transcript}'\
            '?content-type=text/plain&type=cdna')
    return r.text.replace('T', 'U')

beauty = lambda s, x: '\n'.join(splitk(s, x))

if __name__ == '__main__':
    #seq = ''.join([random.choice('AUCG') for k in range(200)])
    with open('mytranscripts.txt', 'r') as f:
        trs = f.read().split('\n')
    for tr in trs:
        print(tr)
        seq = get_sequence(tr)
        db, chosen = p_mat_conf(seq)
        wrt_str = '[SEQUENCE]\n'
        wrt_str += beauty(seq, 80)
        wrt_str += '\n'
        wrt_str += '[db_rndc]\n'
        wrt_str += beauty(db, 80)
        wrt_str += '\n'
        with open(f'results/{tr}/dot_bracket-db_rndc.txt', 'w') as f:
            f.write(wrt_str)
        RNA.svg_rna_plot(seq, db, f'results/{tr}/{tr}-db_rndc.svg')