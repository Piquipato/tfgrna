import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
import random, sys, os
import multiprocessing as mp
from textwrap import wrap as splitk
import time
import tqdm
import RNA

def sph(r, theta, phi): # Generator of Spherical Coordinate Vector
    return (r*np.array([np.sin(theta)*np.cos(phi), \
                        np.sin(theta)*np.sin(phi), \
                        np.cos(theta)]))

def unpack_vector(X):
    x, y, z, xp, yp, zp = X
    X0 = np.array([x, y, z])
    X1 = np.array([xp, yp, zp])
    a = X1 - X0
    return X0, X1, a

def check_sequences(seqA, seqB):
    if len(seqA) < len(seqB):
        null = []
        for i in range(len(seqA), len(seqB)):
            null.append('0')
        seqA += ''.join(null)
    elif len(seqA) > len(seqB):
        null = []
        for i in range(len(seqB), len(seqA)):
            null.append('0')
        seqB += ''.join(null)
    if len(seqA) != len(seqB):
        raise RuntimeError('Sequences not aligned! Check for errors!')
    return seqA, seqB

def check_binds(A, B, nt):
    A0, A1, a = unpack_vector(A)
    B0, B1, b = unpack_vector(B)
    midA = (A0 + A1) / 2
    midB = (B0 + B1) / 2
    dcnt = midB - midA
    if np.dot(dcnt, dcnt) <= (np.dot(a, a) + nt ** 2):
        return True
    return False
    
def binding_links(bdna):
    allowed_binds = [
        'AU',
        'UA',
        'AT',
        'TA',
        'CG',
        'GC',
        'GU',
        'UG'
    ]
    for bind in allowed_binds:
        if bind in bdna:
            return True
    return False

def find_binds(A, seqA, nA, B, seqB, nB, nt):
    binds = []
    A0, A1, a = unpack_vector(A)
    B0, B1, b = unpack_vector(B)
    seqA, seqB = check_sequences(seqA, seqB)
    bps = len(seqA)
    basesA = []
    basesB = []
    for i in np.arange(0, bps):
        basesA.append(A0 + (i / (bps - 1)) * a)
        basesB.append(B0 + (i / (bps - 1)) * b)
    for i in range(len(basesA)):
        for j in range(len(basesB)):
            if np.dot(basesA[i] - basesB[j], basesA[i] - basesB[j]) <= (nt ** 2):
                link = f'b.{nA}.{i}>{nB}.{j}:{seqA[i]}{seqB[j]}'
                if binding_links(link):
                    binds.append(link)
    return binds

"""
def find_binds(data_a, data_b, nt):
    return find_binds(data_a[0], data_a[1], data_a[2], \
                      data_b[0], data_b[1], data_b[2], nt)
"""
                      
def find_seqs(seq, sgs, bps):
    r = len(seq) % bps
    seqs = []
    for i in range(sgs):
        curr = ''
        for j in range(bps):
            try:
                curr += seq[bps*i + j].upper()
            except IndexError:
                curr += '0'
        seqs.append(curr)
    return seqs

def store_list(lst, file='list.txt'):
    with open(file, 'w') as f:
        f.write('\n'.join(lst))

def split_list(lst, splits):
    elements = int(np.floor(len(lst) / splits))
    lstofsplits = []
    for i in range(splits):
        thesplit = []
        for j in range(elements):
            thesplit.append(lst[splits*i + j])
        lstofsplits.append(thesplit)
    return lstofsplits

def loading(current, total):
  perc = (current / total) * 100
  prog = int(np.floor(perc / 5))
  sys.stdout.write('\r')
  sys.stdout.write("[%-20s] %d%% %d / %d" % ('='*prog, perc, current, total))
  sys.stdout.flush()

def get_binds(tpl_ij):
    i, j, vec_i, vec_j, seq_i, seq_j, nt = tpl_ij
    if check_binds(vec_i, vec_j, nt):
        current = find_binds(vec_i, seq_i, i, \
                             vec_j, seq_j, j, nt)
        return current
    return []

def calc_bptbl(vecs, seqs, log=False):
    indeces = []
    coords = []
    nucleotides = []
    positions = []
    k = 0
    for i in range(len(seqs)):
        if log:
            loading(i+1, len(seqs))
        seq = seqs[i]
        X0, X1, x = unpack_vector(vecs[i, :])
        for j in range(len(seq)):
            sumvec = (x * (j / (len(seq) - 1)))
            positions.append(X0 + sumvec)
            indeces.append(k)
            nucleotides.append(seq[j])
            coords.append(f'{i}.{j}')
            k += 1
    data_dict = {'Index': indeces, 'Nucleotide': nucleotides, \
                 'Position': positions, 'Coordinates': coords}
    return pd.DataFrame(data=data_dict)

def calc_fit(arg_tpl):
    i, j, tbl, nt = arg_tpl
    index = list(tbl['Index'])
    pos1 = tbl.loc[index[i], 'Position']
    pos2 = tbl.loc[index[j], 'Position']
    crd1 = tbl.loc[index[i], 'Coordinates']
    crd2 = tbl.loc[index[j], 'Coordinates']
    bp1 = tbl.loc[index[i], 'Nucleotide']
    bp2 = tbl.loc[index[j], 'Nucleotide']
    find = ''
    if binding_links('{bp1}{bp2}'.format(
        bp1=bp1,
        bp2=bp2
    )):
        D = pos1 - pos2
        d = np.sqrt(np.dot(D, D))
        if d < nt:
            what = '{bp1}{bp2}'.format(bp1=bp1, bp2=bp2).upper()
            if what == 'GU' or what == 'UG':
                pre = 'nb.'
            else:
                pre = 'b.'
            find = '{pre}{coord1}>{coord2}:{bp1}{bp2}'.format(
                pre=pre,
                coord1=crd1,
                coord2=crd2,
                bp1=bp1,
                bp2=bp2
            )
    return find

def old_match_binds(tbl, nt):
    print('[INFO] : Entered into old_match_binds w/no loading')
    index = list(tbl['Index'])
    binds = []
    k = 0
    for i in range(len(index) - 1):
        #loading(i+1, len(index) - 1)
        for j in range(i+1, len(index)):
            #loading(k+1, len(index) * (len(index) - 1) / 2)
            if binding_links('{bp1}{bp2}'.format(
                        bp1=tbl.loc[index[i], 'Nucleotide'],
                        bp2=tbl.loc[index[j], 'Nucleotide']
            )):
                D = tbl.loc[index[i], 'Position'] - \
                    tbl.loc[index[j], 'Position']
                d = np.sqrt(np.dot(D, D))
                if d < nt:
                    find = 'b.{coord1}>{coord2}:{bp1}{bp2}'.format(
                        coord1=tbl.loc[index[i], 'Coordinates'],
                        coord2=tbl.loc[index[j], 'Coordinates'],
                        bp1=tbl.loc[index[i], 'Nucleotide'],
                        bp2=tbl.loc[index[j], 'Nucleotide']
                    )
                    binds.append(find)
            k += 1
    return binds

def optimized_match_binds(tbl, bps, nt, log=False):
    if log:
        print('[INFO] : Entered into optimized_match_binds w/no loading')
    index = tbl['Index'].values
    binds = []
    positions = tbl['Position'].values
    coords = tbl['Coordinates'].values
    nucleotides = tbl['Nucleotide'].values

    num_rows = len(index)
    k = 0
    bound = []
    for i in range(num_rows - 1):
        #loading(i+1, num_rows - 1)
        for j in range(i + 1, num_rows):
            #loading(k+1, num_rows * (num_rows - 1) / 2)
            if binding_links(nucleotides[i] + nucleotides[j]) and \
                (not index[i] in bound) and \
                (not index[j] in bound):
                D = positions[i] - positions[j]
                d = np.sqrt(np.dot(D, D))
                if ((d < (2 * nt))) and \
                    (np.abs(index[j] - index[i]) >= np.max([4, bps])):
                    find = 'b.{coord1}>{coord2}:{bp1}{bp2}'.format(
                        coord1=coords[i],
                        coord2=coords[j],
                        bp1=nucleotides[i],
                        bp2=nucleotides[j]
                    )
                    bound.append(index[i])
                    bound.append(index[j])
                    binds.append(find)
            k += 1

    return binds

def match_binds(tbl, nt, multi=True, mtch=False):
    index = list(tbl['Index'])
    states = []
    binds = []
    k = 0
    for i in range(len(index) - 1):
        for j in range(i+1, len(index)):
            loading(k+1, len(index) * (len(index) - 1) / 2)
            states.append((
                i,
                j,
                tbl,
                nt
            ))
            """
                if binding_links('{bp1}{bp2}'.format(
                            bp1=tbl.loc[index[i], 'Nucleotide'],
                            bp2=tbl.loc[index[j], 'Nucleotide']
                )):
                    D = tbl.loc[index[i], 'Position'] - \
                        tbl.loc[index[j], 'Position']
                    d = np.sqrt(np.dot(D, D))
                    if d < nt:
                        find = 'b.{coord1}>{coord2}:{bp1}{bp2}'.format(
                            coord1=tbl.loc[index[i], 'Coordinates'],
                            coord2=tbl.loc[index[j], 'Coordinates'],
                            bp1=tbl.loc[index[i], 'Nucleotide'],
                            bp2=tbl.loc[index[j], 'Nucleotide']
                        )
                        binds.append(find)
            """
            k += 1      
    if not multi:
        binds = []
        k = 1
        for state in states:
            loading(k, len(states))
            find = calc_fit(state)
            if find != '':
                binds.append(find)
            k += 1
    if multi:
        with mp.Pool() as pool:
            """
            binds = list(tqdm.tqdm(pool.imap_unordered(calc_fit, states), \
                                total=len(states)))
            """
            binds = []
            for i, x in enumerate(pool.imap_unordered(calc_fit, states)):
                binds.append(x)
                loading(i+1, len(states))
        """
        with mp.Pool() as pool:
            binds = pool.map(calc_fit, states)
        """
        binds = [it for it in binds if it != '']
        
    return binds

def find_bindseq(seq, tbl, binds):
    bindseq = ['.' for i in range(len(seq))]
    bindlst = []
    index = list(tbl['Index'].values)
    crds = list(tbl['Coordinates'].values)
    bindmat = np.zeros((len(seq), len(seq)))
    for bind in binds:
        crd1 = bind.split('>')[0][2:]
        crd2 = bind.split('>')[1].split(':')[0]
        idx1 = index[crds.index(crd1)]
        idx2 = index[crds.index(crd2)]
        bindmat[crds.index(crd1), crds.index(crd2)] += 1
        bindmat[crds.index(crd2), crds.index(crd1)] += 1
        bindseq[idx1] = '('
        bindseq[idx2] = ')'
        bindlst.append((idx1, idx2))
    return ''.join(bindseq), bindlst, bindmat


def conformate(seq, lp, nt=0.34, processing=False, method='new', mtch=False,
               plotbp=False, plot=True, log=False, path='./out.png'):
    N = len(seq)
    mypath = os.path.dirname(path)
    basename = os.path.basename(path)[:-4]
    if log:
        print(np.ceil(nt * N / (2*lp)), np.ceil((2*lp) / nt))
        print(N, np.ceil(nt * N / (2*lp)) * np.ceil((2*lp) / nt))
    #sgs = int(np.ceil(nt * N / (2*lp)))  number of segments
    bps = int(np.ceil((2*lp) / nt)) # nucleotides per segment
    sgs = int(np.ceil(N / bps)) # number of segments
    #print(len(seq), bps*(sgs-1) + len(seq) % bps)

    sphs = []
    seqs = find_seqs(seq, sgs, bps)
    for k in range(int(sgs)):
        sphs.append(sph(2*lp, np.pi*random.random(), 2*np.pi*random.random()))
    vecs = []
    vecsplot = []
    sumsph = np.zeros_like(sphs[0])
    for k in range(int(sgs)):
        vecs.append(list(np.array([sumsph, sumsph + sphs[k]]).flatten()))
        vecsplot.append(list(np.array([sumsph, sphs[k]]).flatten()))
        sumsph += sphs[k]
    vecs = np.array(vecs, dtype=np.double)
    lim = 1.05*np.max([np.abs(np.min(vecs)), np.abs(np.max(vecs))])
    if log:
        print(len(vecs), sgs)

    binds = []
    if method == 'old':
        if not processing:
            with open('binds.txt', 'a') as f:
                start_time = time.perf_counter()
                for i in range(vecs.shape[0] - 1):
                    if log:
                        loading(i, vecs.shape[0] - 2)
                    for j in range(i+1, vecs.shape[0]):
                        if check_binds(vecs[i, :], vecs[j, :], nt):
                            current = find_binds(vecs[i], seqs[i], i, \
                                                vecs[j], seqs[j], j, nt)
                            for bind in current:
                                f.write(bind + '\n')
                                binds.append(bind)
                finish_time = time.perf_counter()
        else:
            states = []
            for i in range(vecs.shape[0] - 1):
                for j in range(i+1, vecs.shape[0]):
                    states.append((i, j, vecs[i, :], vecs[j, :], \
                                seqs[i], seqs[j], nt))
            start_time = time.perf_counter()
            with mp.Pool() as pool:
                results = pool.map(get_binds, states)
            finish_time = time.perf_counter()
            for result in results:
                if len(result) != 0:
                    for bind in result:
                        binds.append(bind)
    if method == 'new':
        start_time = time.perf_counter()
        tbl = calc_bptbl(vecs, seqs, log=log)
        tbl.to_excel(mypath + f'/{basename}-nucleotides.xlsx', index=True)
        if log:
            sys.stdout.write('\n')
        if mtch:
            binds = match_binds(tbl, nt, multi=processing)
        else:
            binds = optimized_match_binds(tbl, bps, nt, log=log)
        if log:
            sys.stdout.write('\n')
        finish_time = time.perf_counter()
    if log:
        print('Program finished in {secs} seconds with mp: ' \
            '{multi}'.format(secs=(finish_time-start_time), multi=processing))
        print(len(seqs), len(seqs[0]), len(seqs[-1]), len(seq) % bps)
        print(np.sum([len(seqs[i]) for i in range(len(seqs))]), len(seq))
        print(seq[-(len(seq) % bps):])
        print(len(binds) / (2 * len(seq) - 5) * 100)
        print(len(binds) / (np.ceil(N/2) - 1) * 100)
    store_list(binds, file=mypath + f'/{basename}-bindings.txt')
    store_list(seqs, mypath + f'/{basename}-sequences.txt')
    bindseq, bindlst, bindmat = find_bindseq(seq, tbl, binds)
    beauty = lambda s, x: '\n'.join(splitk(s, x))
    bindseqout = '[INFO]: Main sequence of RNA: \n' + beauty(seq, 70) + '\n' + \
        '[INFO]: Binding diagram of RNA: \n' + beauty(bindseq, 70) + '\n'
    with open(mypath + f'/{basename}-bindseq.txt', 'w') as f:
        f.write(bindseqout)

    if plot:
        X, Y, Z, U, V, W = zip(*vecsplot)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        """
        for k in range(len(X)):
            ax.arrow3D(X[k], Y[k], Z[k],
                    U[k], V[k], W[k],
                    arrowstyle="-|>")
        
        """
        ax.quiver(X, Y, Z, U, V, W,
                linewidth=5*nt, edgecolors="#cc0000")
        ax.set_xlim([-lim,lim])
        ax.set_ylim([-lim,lim])
        ax.set_zlim([-lim,lim])
        if method == 'new' and plotbp:
            positions = tbl['Position'].values
            for position in positions:
                ax.scatter(position[0], position[1], position[2], \
                        marker='o', facecolors='#cc0000')
        fig.savefig(path)
        plt.close()
    return binds, bindseq, bindlst, bindmat

def mult_conformate(seq, **kwargs):
    nt = kwargs['nt'] if 'nt' in kwargs else 0.34
    lp = kwargs['lp'] if 'lp' in kwargs else (3/2)*nt
    path = kwargs['path'] if 'path' in kwargs else './out.png'
    plot = kwargs['plot'] if 'plot' in kwargs else False
    iter = kwargs['iter'] if 'iter' in kwargs else 10
    matrix = np.zeros((len(seq), len(seq)))
    for i in range(iter):
        binds, bindseq, bindlst, bindmat = conformate(
            seq, lp, processing=False, nt=nt, plot=plot, path=path
        )
        matrix += bindmat
    matrix /= np.sum(matrix)
    #matrix /= iter
    return matrix


if __name__ == '__main__':
    with open('example.txt', 'r') as f:
        sequence = ''.join(f.read().split('\n'))
    #sequence = 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG' #ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG'
    #sequence = ''.join([random.choice('AUCG') for i in range(random.randint(0, 5000))])
    mp.set_start_method('spawn')
    print(len(sequence))
    nt=0.34
    lp = (3/2)*nt #20.5
    sims = []
    for i in range(10):
        sims.append(conformate(sequence, lp, processing=False, nt=nt, plot=False))