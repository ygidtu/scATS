import numpy as np
from ATS_core import AtsModel, est_density, Parameters, plot_est_vs_real
from ATS_core import Window, WinList, assign_isoform, map_iso_to_gene, build_tree, query, map_to_gene
from csv import reader

def readcsv(filename):
    # read csv file as a list of lists
    with open(filename, 'r') as read_obj:
        # pass the file sobject to reader() to get the reader object
        csv_reader = reader(read_obj, delimiter=' ')
        # Pass reader object to list() to get a list of lists
        list_of_rows = list(csv_reader)
        list_of_rows = [[int(i) for i in l] for l in list_of_rows]
        return list_of_rows

def output_reads_location(filename, list_winlist):
    with open(filename, 'w') as fh:
        for winlist in list_winlist:
            tmplist = ['{}{}{}'.format(win.start, '\t', win.end) for win in winlist]
            fh.write("\t".join(tmplist)+"\n")


######### example 1: ATS mixture #######
def ATS_mixture_example():
    # simulation
    mu_f = 350
    sigma_f = 30
    L = 1500
    n_frag = 1000
    f_len_arr = np.rint(np.random.normal(loc=mu_f, scale=sigma_f, size=n_frag))

    alpha_arr = np.sort(np.array([500, 1000, 1500]))
    n_ats = len(alpha_arr)
    beta_arr = np.random.choice([10, 20, 30], size=n_ats, replace=True)

    # generate weights
    # eg. unif_ws = 0.05  ws = [0.28863634 0.43704733 0.22431633]
    unif_ws = 0.05
    ws = 1 + np.random.uniform(size=n_ats)
    ws = ws / sum(ws)
    ws = (1 - unif_ws) * ws

    boarder = np.rint(n_frag * np.cumsum(ws))  # eg.  [289. 726. 950.]
    seg_st = np.insert(boarder, 0, 0)  # eg.  [  0. 289. 726. 950.]
    seg_en = np.append(boarder, n_frag)  # eg.  [ 289.  726.  950. 1000.]

    label_arr = np.zeros(n_frag, dtype='int')
    st_arr = np.zeros(n_frag, dtype='int')
    en_arr = np.zeros(n_frag, dtype='int')

    for i in range(n_ats + 1):
        tmpinds = np.arange(seg_st[i], seg_en[i], dtype='int')

        tmpn = len(tmpinds)
        label_arr[tmpinds] = i
        # ATS component
        if (i < n_ats):
            st_arr[tmpinds] = np.rint(np.random.normal(loc=alpha_arr[i], scale=beta_arr[i], size=tmpn))
            en_arr[tmpinds] = st_arr[tmpinds] + f_len_arr[tmpinds]
        # uniform component
        else:
            st_arr[tmpinds] = np.random.choice(range(L), size=tmpn, replace=True)
            en_arr[tmpinds] = np.random.choice(range(L), size=tmpn, replace=True)

    atsmix = AtsModel(n_max_ats=7, n_min_ats=1, st_arr=st_arr, gene_length=L,
                 mu_f=mu_f,  # fragment length mean
                 sigma_f=sigma_f,  # fragment length standard deviation

                 # pa site information
                 min_ws=0.05,  # minimum weight of ATS site

                 # inference with fixed parameters
                 fixed_inference_flag=False,

                 debug=False)
    res = atsmix.run()

    # infer with given parameter
    res = atsmix.fixed_inference(res)

    fixed_para = Parameters(title='ground truth', alpha_arr=alpha_arr, beta_arr=beta_arr, L=L)
    res = atsmix.fixed_inference(fixed_para)

    ws = np.append(ws, 0)
    real_para = Parameters(title='ground truth', alpha_arr=alpha_arr, beta_arr=beta_arr, ws=ws, L=L)

    # visualize the estimated density versus real density
    plot_est_vs_real(res, real_para)

    print()
    print(20 * '*' + ' Ground Truth ' + 20 * '*')
    print("real ws: ", ws)
    print("real alpha: ", alpha_arr)
    print("real beta: ", beta_arr)


######### example 2: isoform inference given ATS #######
def isoform_infer_example():
    # prepare data
    filename = 'data/iso.txt'  # isoform information, columns (iso_id, start, end), iso_id=0 is the whole gene
    iso_tbl = readcsv(filename)

    filename = 'data/R1.txt'  # read 1 information, (iso_id, start, end)
    r1_tbl = readcsv(filename)

    filename = 'data/R2.txt'  # read 2 information, (iso_id, start, end)
    r2_tbl = readcsv(filename)

    # break all isoforms into atomic non-overlapping windows
    all_wins = WinList([Window(row[1], row[2]) for row in iso_tbl])
    all_win_on_gene = all_wins.split()
    # print(all_win_on_gene)

    # map all isoforms to the gene (atomic non-overlapping windows)
    iso_inds = list(set([row[0] for row in iso_tbl]))
    n_iso = len(iso_inds)
    iso_wins_list = []
    for i in range(n_iso):
        tmp_wins = WinList([Window(row[1], row[2]) for row in iso_tbl if row[0] == i])
        tmp_wins.sort()
        iso_wins_list.append(map_iso_to_gene(tmp_wins, all_win_on_gene))

    # build the tree
    tree = build_tree(iso_wins_list)

    # test case 1
    testwinlist = WinList([Window(2000, 2032), Window(2572, 3000)])
    res = query(tree, testwinlist, 0)
    assert (res[1])

    # test case 2
    testwin = Window()
    assert (testwin == testwin.shift_start(-100))

    # test case 3
    # query the position of the reads
    n_frag = len(r1_tbl)
    r1st = np.zeros((n_frag, n_iso), dtype=int)
    r1en = np.zeros((n_frag, n_iso), dtype=int)
    r2st = np.zeros((n_frag, n_iso), dtype=int)
    r2en = np.zeros((n_frag, n_iso), dtype=int)

    for i in range(n_frag):
        # print(i)
        e = r1_tbl[i]
        r1_wins = query(tree, WinList([Window(e[1], e[2])]), e[0])
        e = r2_tbl[i]
        r2_wins = query(tree, WinList([Window(e[1], e[2])]), e[0])
        for j, (w1, w2) in enumerate(zip(r1_wins, r2_wins)):
            tmpflag = w1.is_empty() or w2.is_empty() or w1.start < 0 or w2.start < 0
            if tmpflag:
                continue
            r1st[i, j], r1en[i, j], r2st[i, j], r2en[i, j] = w1.start, w1.end, w2.start, w2.end

    R1ST = np.genfromtxt('data/R1ST.txt', delimiter=' ')
    R1ST[R1ST == -1] = 0
    assert np.all(R1ST == r1st[:, 1:])

    R2ST = np.genfromtxt('data/R2ST.txt', delimiter=' ')
    R2ST[R2ST == -1] = 0
    assert np.all(R2ST == r2st[:, 1:])

    R1EN = np.genfromtxt('data/R1EN.txt', delimiter=' ')
    assert np.all(R1EN == r1en[:, 1:])

    R2EN = np.genfromtxt('data/R2EN.txt', delimiter=' ')
    assert np.all(R2EN == r2en[:, 1:])

    # test case 4
    # estimate the proportion of different isoforms given an ATS
    ats_pos = 1400  # 1400, 2500, 4000, 1000
    mu_frag = 350
    sd_frag = 50
    # r1_wins_list = [WinList([Window(e[1], e[2]), ]) for e in r1_tbl if e[0] == 1]
    # r2_wins_list = [WinList([Window(e[1], e[2]), ]) for e in r2_tbl if e[0] == 1]
    # read_labels = [1 for e in r1_tbl if e[0]==1]

    r1_wins_list = [WinList([Window(e[1], e[2]), ]) for e in r1_tbl]
    r2_wins_list = [WinList([Window(e[1], e[2]), ]) for e in r2_tbl]
    read_labels = [e[0] for e in r1_tbl]

    n_iso = len(iso_inds)
    iso_wins_list = []
    for i in range(n_iso):
        tmp_wins = WinList([Window(row[1], row[2]) for row in iso_tbl if row[0] == i])
        tmp_wins.sort()
        iso_wins_list.append(tmp_wins)

    res = assign_isoform(ats_pos, iso_wins_list, r1_wins_list, r2_wins_list, read_labels, mu_frag, sd_frag)
    res = np.around(res[1:], decimals=3)  # remove weights for the first exon
    if all(res < 1e-8):
        print(f'ats_pos={ats_pos} potentially support a novel isoform')
    print(f'iso_weights={res}')

    # test case 5
    # map a window to gene
    qwinlist = r1_wins_list[0]
    res_winlist = map_to_gene(tree, qwinlist,  read_labels[0])
    for win in res_winlist:
        print(win)

    r1_on_gene = []
    for i,  win in enumerate(r1_wins_list):
        tmp_winlist = map_to_gene(tree, win, read_labels[i])
        r1_on_gene.append(tmp_winlist)
    # output_reads_location('data/r1_on_gene.txt', r1_on_gene)

    r2_on_gene = []
    for i, win in enumerate(r2_wins_list):
        tmp_winlist = map_to_gene(tree, win, read_labels[i])
        r2_on_gene.append(tmp_winlist)
    # output_reads_location('data/r2_on_gene.txt', r2_on_gene)

    # check if the we can map them back to the correct isoforms
    # query the position of the reads
    r1st_new = np.zeros((n_frag, n_iso), dtype=int)
    r1en_new = np.zeros((n_frag, n_iso), dtype=int)
    r2st_new = np.zeros((n_frag, n_iso), dtype=int)
    r2en_new = np.zeros((n_frag, n_iso), dtype=int)

    for i in range(n_frag):
        # print(i)
        r1_wins = query(tree, r1_on_gene[i], 0)
        r2_wins = query(tree, r2_on_gene[i], 0)
        for j, (w1, w2) in enumerate(zip(r1_wins, r2_wins)):
            tmpflag = w1.is_empty() or w2.is_empty() or w1.start < 0 or w2.start < 0
            if tmpflag:
                continue
            r1st_new[i, j], r1en_new[i, j], r2st_new[i, j], r2en_new[i, j] = w1.start, w1.end, w2.start, w2.end
    assert np.all(r1st_new == r1st)
    assert np.all(r1en_new == r1en)

    assert np.all(r2st_new == r2st)
    assert np.all(r2en_new == r2en)




if __name__ == '__main__':
    #ATS_mixture_example()
    isoform_infer_example()
