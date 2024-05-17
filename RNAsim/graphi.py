import RNA


if __name__ == '__main__':
    with open('mytranscripts.txt') as g:
        trs = g.read().split('\n')

    for tr in trs:
        db_path = f'results/{tr}/dot_bracket.txt'
        db_path_rndc = f'results/{tr}/dot_bracket-db_rndc.txt'
        seq_path = f'results/{tr}/basesequence.txt'
        print(f'[{tr}] [{db_path}][{seq_path}]:')
        with open(db_path, 'r') as f:
            db = [it for it in f.read().split('\n') if 'db_one' in it][0].split('\t')[0]
        with open(seq_path, 'r') as f:
            seq = f.read()
        print(seq)
        print(db)
        linked_ncl = len(db) - len([it for it in list(db) if it == '.'])
        print(linked_ncl, -linked_ncl/2)
        RNA.svg_rna_plot(seq, db, f'results/{tr}/{tr}-db_one_conf.svg')
        with open(db_path_rndc, 'r') as f:
            db = f.read()
        linked_ncl = len(seq) - len([it for it in list(db) if it == '.'])