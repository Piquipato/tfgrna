import numpy as np

def calc_link(s):
    n = len(s)
    link = np.zeros((n, n))
    queue = []
    for i in range(len(s)):
        if s[i] == '(':
            queue.append(i)
        if s[i] == ')':
            j = queue.pop()
            link[i, j] = 1
            link[j, i] = 1
    return link

def calc_Zij(link, i, j, T=37, R=8.314):
    if np.abs(i-j) < 4:
        with open('zij.txt', 'a') as f:
            f.write(f'{i},{j},1\n')
        return 1
    db = {}
    with open('zij.txt', 'r') as f:
            calcs = f.read().split('\n')
            for calc in calcs:
                if calc != '':
                    (k,l,val) = calc.split(',')
                    db[(k, l)] = float(val)
                    db[(k, l)] = float(val)
    if (i, j) in db.keys():
        return db[(i, j)]
    else:
        sum = calc_Zij(link, i, j-1, T, R)
        for k in range(i, j-4):
            sum += calc_Zij(link, i, k-1, T, R) * calc_Zij(link, k+1, j-1, T, R) * np.exp(-link[k, j]/(R*T))
        with open('zij.txt', 'a') as f:
            f.write(f'{i},{j},{sum}\n')
        return sum

def calc_Z(s, T=37, R=8.314):
    link = calc_link(s)
    n = len(s)
    Z = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            val = calc_Zij(link, i, j, T, R)
            Z[i, j] = val
            Z[j, i] = val
    return Z

def calc_pij(s, Z, T=37, R=8.314):
    link = calc_link(s)
    n, m = np.shape(Z)
    p = np.zeros(np.shape(Z))
    for i in range(n):
        for j in range(m):
            try:
                p[i, j] = (Z[0, i-1] * Z[j+1, n-1] * np.exp(- link[i, j] / (R*T))) / Z[0, n-1]
            except:
                p[i, j] = 0
    return p

def calc_nussinov(s, T=37, R=8.314):
    n = len(s)
    Z = calc_Z(s, T, R)
    p = calc_pij(s, Z, T, R)
    return -np.log(Z[0, n-1])

if __name__ == '__main__':
    import random
    import RNA
    seq = ''.join([random.choice('AUCG') for i in range(200)])
    (ss, mfe) = RNA.fold_compound(seq).mfe()
    Z = calc_Zij(calc_link(ss), 0, len(ss) - 1)
    print(-np.log(Z))