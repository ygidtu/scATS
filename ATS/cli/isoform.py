#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2021.05.06 by Zhang
"""

@click.command()
@click.option(
    "--utr",
    type=click.Path(exists = True),
    required=True,
    help=""" The path to utr file, bed format. """
)
@click.option(
    "-o", "--output", 
    type=click.Path(),
    required=True,
    help=""" The path to output file. """
)
@click.option(
    "--n-max-ats",
    type=click.IntRange(1, math.inf),
    default = 5,
    help=""" The maximum number of ATSs in same UTR. """
)
@click.option(
    "--n-min-ats",
    type=click.IntRange(1, math.inf),
    default = 1,
    help=""" The minimum number of ATSs in same UTR. """
)
@click.option(
    "--utr-length",
    type=int,
    default = 1000,
    help=""" The length of UTR. """
)
@click.option(
    "--utr-length",
    type=int,
    default = 1000,
    help=""" The estimate length of gene. """
)
@click.option(
    "--mu-f",
    type=int,
    default = 300,
    help=""" The mean of fragment length. """
)   
@click.option(
    "--sigma-f",
    type=int,
    default = 50,
    help=""" The standard deviation of fragment length. """
)
@click.option(
    "--min-ws",
    type=float,
    default = 0.01,
    help=""" The minimum weight of ATSs. """
)
@click.option(
    "--max-unif-ws",
    type=float,
    default = 0.1,
    help=""" The maximum weight of uniform component. """
)
@click.option(
    "--max-beta",
    type=int,
    default = 50,
    help=""" The maximum std for ATSs. """
)
@click.option(
    "--fixed-inference",
    is_flag=True,
    default = True,
    type=click.BOOL,
    help=""" Inference with fixed parameters. """
)
@click.option(
    "-d", "--debug",
    is_flag=True,
    type=click.BOOL,
    help=""" Enable debug mode to get more debugging information. """
)
@click.option(
    "-p", "--processes",
    type=click.IntRange(1,cpu_count()),
    default = 1,
    help=""" How many cpu to use. """
)
@click.argument("bams", nargs = -1, type=click.Path(exists=True), required=True)
def inference(
    utr: str,
    output: str,
    n_max_ats: int, 
    n_min_ats: int,
    utr_length: int,
    mu_f: int,
    sigma_f: int,
    min_ws: float,
    max_unif_ws: float,
    max_beta: int,
    fixed_inference: bool,
    processes: int,
    debug: bool, 
    bams: List[str],
):
    u"""
    Inference
    \f

    :param debug: enable debug mode
    """

    init_logger("DEBUG" if debug else "INFO")

    params = ATSParams(
        utr = utr,
        bam = bams,
        n_max_ats = n_max_ats, 
        n_min_ats = n_min_ats,
        utr_length = utr_length,
        mu_f = mu_f,
        sigma_f = sigma_f,
        min_ws = min_ws,
        max_unif_ws = max_unif_ws,
        max_beta = max_beta,
        fixed_inference_flag = fixed_inference
    )

    # prepare data
    filename = 'data/iso.txt'  # isoform information, columns (iso_id, start, end), iso_id=0 is the whole gene
    iso_tbl = readcsv(filename)

    # break all isoforms into atomic non-overlapping windows
    all_wins = WinList([Window(row[1], row[2]) for row in iso_tbl])
    all_win_on_gene = all_wins.split()
    # print(all_win_on_gene)

    # map all isoforms to the gene (atomic non-overlapping windows)
    iso_inds = list(set([row[0] for row in iso_tbl]))
    n_iso = len(iso_inds)
    iso_wins_on_gene_list = []
    iso_wins_list = []
    for i in range(n_iso):
        tmp_wins = WinList([Window(row[1], row[2]) for row in iso_tbl if row[0] == i])
        tmp_wins.sort()
        iso_wins_list.append(tmp_wins)
        iso_wins_on_gene_list.append(map_iso_to_gene(tmp_wins, all_win_on_gene))

    # build the tree
    tree = build_tree(iso_wins_on_gene_list)

    # get the positions of r1 and r2 on gene
    filename = 'data/r1_on_gene.txt'  # read 1 information, (start1, end1), (start2, end2), ...
    r1_list = read_location(filename)  # list of WinList

    filename = 'data/r2_on_gene.txt'  # read 2 information, (start1, end1), (start2, end2), ...
    r2_list = read_location(filename)  # list of WinList

    n_frag = len(r1_list)
    assert n_frag == len(r2_list)

    frag_label = [0 for _ in range(n_frag)]  # isoform index of fragments, 0 is the whole gene


    # # query the position of r1 and r2 on different isoforms
    # r1st = np.zeros((n_frag, n_iso), dtype=int)
    # r1en = np.zeros((n_frag, n_iso), dtype=int)
    # r2st = np.zeros((n_frag, n_iso), dtype=int)
    # r2en = np.zeros((n_frag, n_iso), dtype=int)
    # st_arr = np.zeros(n_frag, dtype=int)

    # for i in range(n_frag):
    #     # print(i)
    #     st_arr[i] = r1_list[i][0].start
    #     r1_wins = query(tree, r1_list[i], 0)
    #     r2_wins = query(tree, r2_list[i], 0)
    #     for j, (w1, w2) in enumerate(zip(r1_wins, r2_wins)):
    #         tmpflag = w1.is_empty() or w2.is_empty() or w1.start < 0 or w2.start < 0
    #         if tmpflag:
    #             continue
    #         r1st[i, j], r1en[i, j], r2st[i, j], r2en[i, j] = w1.start, w1.end, w2.start, w2.end

    # prepare results
    data = PreprocData()
    # data.st_arr = st_arr
    # data.r1st = r1st
    # data.r1en = r1en
    # data.r2st = r2st
    # data.r2en = r2en

    data.st_arr = np.array([winlist[0].start for winlist in r1_list])

    data.iso_wins_on_gene_list = iso_wins_on_gene_list
    data.iso_wins_list = iso_wins_list
    data.tree = tree
    data.n_frag = n_frag
    data.n_iso = n_iso
    data.gene_length = all_win_on_gene[-1].end

    data.r1_list = r1_list
    data.r2_list = r2_list
    data.frag_label = frag_label

    data.n_max_ats = 5
    data.n_min_ats = 1
    data.mu_f = 500
    data.sigma_f = 30
    data.min_ws = 0.05
    data.max_unif_ws = 0.1

    return data