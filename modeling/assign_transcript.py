#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from __future__ import annotations

import os
import re

from functools import lru_cache, reduce
from itertools import chain
from typing import List

import filetype
import numpy as np
import pysam

from loguru import logger
from rich import print


class Window:
    DEL_VAL = 1 << 31

    def __init__(self, start: int = DEL_VAL, end: int = DEL_VAL) -> Window:
        """
        :type start: int
        :type end: int
        """
        self.start = start
        self.end = end

    def __key(self):
        return self.start, self.end

    def __hash__(self):
        return hash(self.__key())

    def is_empty(self):
        return self.start == self.end

    def __bool__(self):
        return self.start != Window.DEL_VAL

    def __len__(self):
        return self.end - self.start

    def __str__(self):
        return f'{self.start}  {self.end}  {len(self)}'

    #################### window vs window -> bool #####################
    def __lshift__(self, other):
        return self.end <= other.start

    def __lt__(self, other):
        return self.start < other.start < self.end and other.start < self.end < other.end

    def __le__(self, other):
        return other.start <= self.start and self.end <= other.end

    def __eq__(self, other):
        if isinstance(other, Window):
            return self.start == other.start and self.end == other.end
        else:
            return False

    def __ge__(self, other):
        return self.start <= other.start and other.end <= self.end

    def __gt__(self, other):
        return other.start < self.start < other.end < self.end

    def __rshift__(self, other):
        return other.end <= self.start

    # overlap
    def __ne__(self, other):
        return not (other.end <= self.start or self.end <= other.start)

    # adjacent
    def adj(self, other):
        return self.end == other.start or self.start == other.end

    #################### window vs point operation #####################
    # if a point falls within the window
    def __contains__(self, point):
        if self.is_empty():
            return False
        return self.start <= point < self.end

    def shift_start(self, start):
        if self.start == Window.DEL_VAL:
            return self
        else:
            return Window(start, start + self.end - self.start)


class WinList(list):
    def __init__(self, *args, is_sorted=False):
        super().__init__(*args)
        self.is_sorted = is_sorted

    def __str__(self):
        return "\n".join(['start  end  length'] + [str(w) for w in self.__iter__()])

    def append(self, __object: Window) -> None:
        super().append(__object)
        self.is_sorted = False

    def rm_empty_win(self) -> WinList:
        return WinList([w for w in self if w], is_sorted=self.is_sorted)

    def sort(self) -> None:
        if not self.is_sorted:
            super().sort(key=lambda win: win.start)
            self.is_sorted = True

    def rm_duplicate(self) -> WinList:
        res = WinList(set(self))
        res.sort()
        return res

    # merge overlapping windows
    def merge(self) -> WinList:
        winlist = self.rm_empty_win()
        assert len(winlist) > 0
        if not winlist.is_sorted:
            winlist.sort()
        res_list = WinList([winlist[0]])
        for win in winlist:
            curr_win = res_list[-1]
            if curr_win.start <= win.start <= curr_win.end:
                res_list[-1] = Window(curr_win.start, win.end)
            else:
                res_list.append(win)
        res_list.is_sorted = True
        return res_list


    # split the winList into tiny non-overlapping intervals
    def split(self):
        winlist = self.rm_empty_win()
        if len(winlist) == 0:
            return winlist
        if not winlist.is_sorted:
            winlist.sort()
        boarder = set()
        for win in winlist:
            if not win:
                continue
            boarder.add(win.start)
            boarder.add(win.end)
        boarder_arr = sorted(list(boarder))
        winlist = [Window(i, j) for i, j in zip(boarder_arr, boarder_arr[1:])]
        return WinList(winlist, is_sorted=True)

    # return the left most position of the windows list
    def get_left(self):
        for win in self:
            if win.start != Window.DEL_VAL:
                return win.start
        return Window.DEL_VAL

    # return the right most position of the windows list
    def get_right(self):
        for win in reversed(self):
            if win.end != Window.DEL_VAL:
                return win.end
        return Window.DEL_VAL

    def get_range(self):
        return Window(self.get_left(), self.get_right())


class TreeNode:
    def __init__(self, left=None, right=None, winlist: WinList = WinList()):
        """
        :type left: TreeNode
        :type right: TreeNode
        :type winlist: a list of windows, each
        """
        self.left = left
        self.right = right
        self.winlist = winlist  # ith window corresponds to range of this tree node


def map_iso_to_gene(iso_list: WinList, gene_list: WinList) -> WinList:
    """
    :param iso_list: exon windows of the isoform on gene
    :param gene_list: smallest windows across isoforms by breaking the gene
    :return: a WinList (same size as gene_list), with indexes relative to the isoform (local index)
    """
    assert iso_list.is_sorted
    assert gene_list.is_sorted
    res_win_list = WinList()
    curr_exon_j = 0
    local_iso_index = 0
    for win in gene_list:
        if curr_exon_j >= len(iso_list):
            res_win_list.append(Window())
            continue
        if win <= iso_list[curr_exon_j]:
            res_win_list.append(win.shift_start(local_iso_index))
            local_iso_index += len(win)
            if iso_list[curr_exon_j].end == win.end:
                curr_exon_j += 1
        else:
            res_win_list.append(Window())
    return res_win_list


def build_tree(isowins_list: List[WinList]):
    """
    :param isowins_list: each element is an aligned WinList of the same window across different isoforms
    :return: a tree
    """

    def dfs(i_left, i_right):
        root = TreeNode()
        if i_left == i_right - 1:
            root.winlist = alignwins_list[i_left]
            root.left = None
            root.right = None
        else:
            mid = (i_left + i_right) // 2
            tmplist = [WinList(isowins[i_left:i_right]).get_range() for isowins in isowins_list]
            root.winlist = WinList(tmplist)
            root.left = dfs(i_left, mid)
            root.right = dfs(mid, i_right)
        return root

    # each element is an aligned WinList of the same window across different isoforms
    alignwins_list = [WinList(wins) for wins in zip(*isowins_list)]
    return dfs(0, len(alignwins_list))


def concatenate_winlist(winlist1: WinList, winlist2: WinList) -> WinList:
    res_list = [concatenate_windows(w1, w2) for w1, w2 in zip(winlist1, winlist2)]
    return WinList(res_list)


def concatenate_windows(win_left: Window, win_right: Window) -> Window:
    if win_left and win_right and win_left.end == win_right.start:
        return Window(win_left.start, win_right.end)
    else:
        return Window()


def query(root: TreeNode, query_winlist: WinList, iso_index: int) -> WinList:
    """
    map the query_winlist to isoform, get the relative index of the query window on all isoforms
    :param root: root of the tree
    :param query_winlist: winlist of a read mapped to the "iso_index"th isoform
    :param iso_index: index of given isoform, used for locating the window
    :return: mapped windows on all isoforms
    """

    @lru_cache(maxsize=1024)
    def query_window(root: TreeNode, qwin: Window, i_iso: int) -> WinList:
        if root is None or root.winlist[i_iso].is_empty():
            return none_winlist
        if not qwin <= root.winlist[i_iso]:
            return none_winlist
        if root.left is None and root.right is None:  # leaf case
            local_st = qwin.start - root.winlist[i_iso].start
            local_win = Window(local_st, local_st + len(qwin))
            res_winlist = WinList()
            for win in root.winlist:
                if not win:
                    res_winlist.append(Window())
                else:
                    res_winlist.append(local_win.shift_start(win.start + local_st))
            return res_winlist
        # root is inner node, one child might be intron
        if root.left.winlist[i_iso].end != Window.DEL_VAL:
            mid = root.left.winlist[i_iso].end
        else:
            mid = root.right.winlist[i_iso].start

        if qwin.end <= mid:
            return query_window(root.left, qwin, i_iso)
        elif qwin.start >= mid:
            return query_window(root.right, qwin, i_iso)
        else:
            qwin_left = Window(qwin.start, mid)
            qwin_right = Window(mid, qwin.end)
            left_winlist = query_window(root.left, qwin_left, i_iso)
            right_winlist = query_window(root.right, qwin_right, i_iso)
            # print(left_winlist)
            return WinList([concatenate_windows(lw, rw) for lw, rw in zip(left_winlist, right_winlist)])

    n_isoform = len(root.winlist)
    none_winlist = WinList([Window() for _ in range(n_isoform)])

    res_list = [query_window(root, win, iso_index) for win in query_winlist]
    res_winlist = reduce(concatenate_winlist, res_list)

    return res_winlist


# modify exons of isoform such that ATS falls into an exon
def proc_isoform(ats_pos: int, iso_wins: WinList) -> WinList:
    """
    :param iso_wins_list: list of isoforms, 0th element is gene range
    :return: a list of isoforms with exons modified
    """
    # step 1: check if [max(0, ats-50), ats+50] falls in exon
    # step 1.1: if yes => no special treatment
    # step 1.2: if no  => merge [max(0, ats-50), min(ats+200, exon_start)] with the nearest first exon
    ats_win = Window(ats_pos, ats_pos+1)
    i = 0
    while iso_wins[i] << ats_win:
        i += 1
    first_exon_ind = i # first downstream exon

    ats_win = Window(max(0, ats_pos-50), ats_pos+50)
    if ats_win <= iso_wins[first_exon_ind]:
        return iso_wins

    iso_wins.append(Window(max(0, ats_pos-50), min(ats_pos+200, iso_wins[first_exon_ind].start)))
    iso_wins.is_sorted = False
    return iso_wins.merge()   # sorted wins


def get_read_relative_pos(iso_wins_list: List[WinList], r1_wins_list: List[WinList], r2_wins_list: List[WinList], read_labels):
    # break all isoforms into atomic non-overlapping windows
    all_wins = WinList([win for iso_wins in iso_wins_list for win in iso_wins])
    all_win_on_gene = all_wins.split()

    # map all isoforms to the gene (atomic non-overlapping windows)
    atom_isowins_list = []
    for iso_wins in iso_wins_list:
        atom_isowins_list.append(map_iso_to_gene(iso_wins, all_win_on_gene))

    # build the tree
    tree = build_tree(atom_isowins_list)

    # get relative position of R1 and R2 on the gene
    n_frag = len(read_labels)
    n_iso = len(iso_wins_list) - 1
    read_st_mat = np.zeros((n_frag, n_iso), dtype=int)
    read_en_mat = np.zeros((n_frag, n_iso), dtype=int)

    for i in range(n_frag):
        r1_wins = r1_wins_list[i]
        r2_wins = r2_wins_list[i]
        res1 = query(tree, r1_wins, read_labels[i])
        res2 = query(tree, r2_wins, read_labels[i])
        for j in range(n_iso):
            if res1[j + 1] and res2[j + 1]:
                read_st_mat[i, j] = res1[j + 1].start
                read_en_mat[i, j] = res2[j + 1].end

    return read_st_mat, read_en_mat


def assign_isoform(ats_pos: int, iso_wins_list: List[WinList], r1_wins_list: List[WinList], r2_wins_list: List[WinList],
                   read_labels, mu_frag, sd_frag):
    # extend ATS and merge with exons
    # 0th element of iso_wins_list is the whole gene window i.e. [Window(0,gene_length),]
    iso_wins_list = [proc_isoform(ats_pos, iso_wins) for iso_wins in iso_wins_list]

    # get relative positions of the reads on each isoform
    read_st_mat, read_en_mat = get_read_relative_pos(iso_wins_list, r1_wins_list, r2_wins_list, read_labels)

    # assign each pair-end read to different isoforms
    n_frag, n_iso = read_st_mat.shape
    read_len_mat = read_en_mat - read_st_mat
    post_prob_mat = cal_post_prob(mu_frag, sd_frag, read_len_mat)


    #TODO
    # 同一个isoform，短片段和长片段的比例，是不是双峰分布，需要前面的mixture代码

    iso_post_prob = np.sum(post_prob_mat, axis=0)
    map_rate = np.round(np.sum(iso_post_prob))/n_frag  # 至少map到一个isoform上的fragment比例
    iso_post_prob = iso_post_prob/np.sum(iso_post_prob) # 对应于该ATS的各isoform比例

    return map_rate, iso_post_prob


# calculate the posterior probability of each DNA fragment on different isoforms
def cal_post_prob(mu_frag_size: int, sd_frag_size: int, len_iso_mat: np.ndarray):
    """
    :param mu_frag_size: mean of fragment size
    :param sd_frag_size: std of fragment size
    :param len_iso_mat: n_frag x n_iso numpy matrix, each element is the length of DNA fragment on each isoform
    :return: n_frag x n_iso numpy matrix, each element is the posterior of a given fragment on an isoform
    """

    def norm_pdf(x, mean, sd):
        prob_density = (1 / (np.sqrt(2 * np.pi) * sd) ) * np.exp(-0.5 * ((x - mean) / sd) ** 2)
        return prob_density

    max_arr = np.max(len_iso_mat, axis=1)
    val_inds = max_arr > 0

    prob_mat = np.zeros(len_iso_mat.shape, dtype="float")
    prob_mat[val_inds, :] = norm_pdf(len_iso_mat[val_inds, :], mu_frag_size, sd_frag_size)
    prob_mat[len_iso_mat==0] = 0

    prob_mat[val_inds, :] = prob_mat[val_inds, :]/np.sum(prob_mat[val_inds, :], 1)[:, np.newaxis]

    return prob_mat


class BED:

    def __init__(self, chrom: str, start: int, end: int, strand: str, name: str = "", parent: str = ""):
        self.chromosome = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.name = name
        self.parent = parent

    def __str__(self):
        res = f"{self.chromosome}:{self.start}-{self.end}:{self.strand}"
        if self.name:
            res = f"{res}\t{self.name}"

        if self.parent:
            res = f"{res}\t{self.parent}"
        return res

    def __hash__(self):
        return hash((self.chromosome, self.start, self.end, self.strand, self.name, self.parent))


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


def load_ats(path: str):
    beds = []
    header = True
    with open(path) as r:
        for line in r:
            if header:
                header = False
                continue

            line = line.strip().split("\t")
            
            if line[-1] == "NA" or "inf" in line[-1].lower():
                continue

            site = re.split(r"[:-]", line[0])
            strand = site[-1]

            if strand != "+":
                strand = "-"
            chrom, start_pos, end_pos = site[0], int(site[1]), int(site[2])
            alpha = line[4].split(",")

            if len(alpha) > 1:
                try:
                    for x in alpha:
                        if x != "":
                            s = int(float(x))
                            
                            if strand == "+":
                                beds.append(BED(chrom, s, s + 1, strand, line[0], str(len(beds) + 1))) 
                            else:
                                beds.append(BED(chrom, s - 1, s, strand, line[0], str(len(beds) + 1)))
                except Exception as e:
                    print(e)
                    pass

    return beds


def is_gtf(infile):
    u"""
    check if input file is gtf
    :param infile: path to input file
    :return:
    """
    if infile is None:
        return False

    is_gtf = 0
    try:
        if filetype.guess_mime(infile) == "application/gzip":
            is_gtf += 10
            r = gzip.open(infile, "rt")
        else:
            r = open(infile)
    except TypeError as err:
        logger.error("failed to open %s", infile)
        exit(err)

    for line in r:
        if line.startswith("#"):
            continue

        lines = re.split(r"\s+", line)

        if len(lines) < 8:
            break

        if re.search(
            r"([\w-]+ \"[\w.\s\-%,:]+\";? ?)+",
            " ".join(lines[8:])
        ):
            is_gtf += 1

        break

    r.close()

    return is_gtf


def index_gtf(input_gtf, sort_gtf=True, retry=0):
    u"""
    Created by ygidtu
    Extract only exon tags and keep it clean
    :param input_gtf: path to input gtf file
    :param sort_gtf: Boolean value, whether to sort gtf file first
    :param retry: only try to sort gtf once
    :return path to compressed and indexed bgzipped gtf file
    """
    if input_gtf is None:
        return None
    gtf = is_gtf(input_gtf)

    if gtf % 10 != 1:
        raise ValueError("gtf file required, %s seems not a valid gtf file" % input_gtf)

    index = False
    if gtf // 10 > 0:
        output_gtf = input_gtf
    else:
        output_gtf = input_gtf + ".gz"
    if not os.path.exists(output_gtf) or not os.path.exists(output_gtf + ".tbi"):
        index = True

    elif os.path.getctime(output_gtf) < os.path.getctime(output_gtf) or \
            os.path.getctime(output_gtf) < os.path.getctime(output_gtf):
        index = True

    # 2018.12.21 used to handle gtf not sorted error
    if sort_gtf and retry > 1:
        raise OSError("Create index for %s failed, and trying to sort it failed too" % input_gtf)
    elif sort_gtf:
        data = []

        logger.info("Sorting %s" % input_gtf)

        old_input_gtf = input_gtf
        input_gtf = re.sub("\.gtf$", "", input_gtf) + ".sorted.gtf"

        output_gtf = input_gtf + ".gz"

        if os.path.exists(input_gtf) and os.path.exists(output_gtf):
            return output_gtf

        try:
            w = open(input_gtf, "w+")
        except IOError as err:
            w = open("/tmp/sorted.gtf")

        with open(old_input_gtf) as r:
            for line in tqdm(r):
                if line.startswith("#"):
                    w.write(line)
                    continue

                lines = line.split()

                if len(lines) < 1:
                    continue

                data.append(
                    GenomicLoci(
                        chromosome=lines[0],
                        start=lines[3],
                        end=lines[4],
                        strand=lines[6],
                        gtf_line=line
                    )
                )

        for i in sorted(data):
            w.write(i.gtf_line)

        w.close()

    if index:
        logger.info("Create index for %s", input_gtf)
        try:
            pysam.tabix_index(
                input_gtf,
                preset="gff",
                force=True,
                keep_original=True
            )
        except OSError as err:

            if re.search("could not open", str(err)):
                raise err

            logger.error(err)
            logger.error("Guess gtf needs to be sorted")
            return index_gtf(input_gtf=input_gtf, sort_gtf=True, retry=retry + 1)

    return output_gtf


def read_transcripts(gtf_file, region:BED):
    u"""
    Read transcripts from tabix indexed gtf files
    The original function check if the junctions corresponding to any exists exons, I disable this here
    :param gtf_file: path to bgzip gtf files (with tabix index), only ordered exons in this gtf file
    :param region: splice region
    :param retry: if the gtf chromosome and input chromosome does not match. eg: chr9:1-100:+ <-> 9:1-100:+
    :param genome: path to genome fasta file
    :return: SpliceRegion
    """

    res = {}

    if gtf_file is None:
        return res

    if not os.path.exists(gtf_file):
        raise FileNotFoundError("%s not found" % gtf_file)
    print(region)
    try:
        logger.info("Reading from %s" % gtf_file)

        with pysam.Tabixfile(gtf_file, 'r') as gtf_tabix:
            relevant_exons_iterator = gtf_tabix.fetch(
                region.chromosome,
                region.start - 1,
                region.end + 1,
                parser=pysam.asGTF()
            )

            # min_exon_start, max_exon_end, exons_list = float("inf"), float("-inf"),  []
            transcripts = []
            genes = set()
            for line in relevant_exons_iterator:
                try:
                    if line.feature == "transcript" and line.strand == region.strand:
                        transcripts.append(BED(
                            line.contig,
                            line.start,
                            line.end,
                            line.strand,
                            line.transcript_id,
                            line.gene_id
                        ))

                        genes.add(line.gene_id)
                except IndexError as err:
                    logger.error(err)
            
            if len(genes) > 1:
                return res

            for t in transcripts:
                exons = []
                for line in gtf_tabix.fetch(t.chromosome, t.start - 1, t.end + 1, parser=pysam.asGTF()):
                    try:
                        if line.feature == "exon" and line.transcript_id == t.name:
                            exons.append(BED(
                                line.contig,
                                line.start,
                                line.end,
                                line.strand,
                                line.exon_id,
                                line.transcript_id
                            ))

                    except IndexError as err:
                        logger.error(err)
                res[t] = exons

    except ValueError as err:
        logger.warn(err)
    return res


def extract_from_cigar_string(record: pysam.AlignedSegment, mode:str=""):
    u"""
    extract junctions from cigar string

    M	BAM_CMATCH	0
    I	BAM_CINS	1
    D	BAM_CDEL	2
    N	BAM_CREF_SKIP	3
    S	BAM_CSOFT_CLIP	4
    H	BAM_CHARD_CLIP	5
    P	BAM_CPAD	6
    =	BAM_CEQUAL	7
    X	BAM_CDIFF	8
    B	BAM_CBACK	9

    IGV only skip S and I

    :param record: SAMSegment
    :return:
    """
    pos = record.reference_start + 1
    pos_list = []

    skipped_code = (1, 2, 4, 5)
    
    modes = {
        "star": (2, 5, 6),
        "igv": (1, 4),
        "none": ()
    }

    skipped_code = modes.get(mode, skipped_code)

    for i, j in record.cigartuples:

        if i not in skipped_code:
            pos += j

        if i == 3:
            pos_list.append(pos - j)
            pos_list.append(pos - 1)
    
    return pos_list


def assign_reads(reads, transcripts:Dict):
    res = []

    introns = extract_from_cigar_string(reads)

    if introns:
        for key, exons in transcripts.items():
            for i in range(0, len(introns), 2):
                for j in range(0, len(exons), 2):
                    if abs(introns[i] - exons[j].end) < 3 and abs(introns[i+1] - exons[j+1].start) < 3:
                        res.append(key)
                        continue

                    if exons[j].start > introns[i] + 3:
                        break
    else:

        for key, exons in transcripts.items():
            for e in exons:
                if e.start <= reads.reference_start and e.end >= reads.reference_end:
                    res.append(key)
    return res


def load_reads(bam_file, region:BED, transcripts):
    res = {
        "R1": [],
        "R2": [],
        "label": set()
    }

    with pysam.AlignmentFile(bam_file) as r:
        for rec in r.fetch(region.chromosome, region.start, region.end):

            if rec.is_unmapped or rec.is_qcfail or rec.mate_is_unmapped:
                continue

            if rec.is_read1:
                try:
                    mate = r.mate(rec)
                    r1_label = assign_reads(rec, transcripts)
                    r2_label = assign_reads(mate, transcripts)

                    for label in set(r1_label) & set(r2_label):
                        res["R1"].append([label, rec.reference_start, rec.reference_end])
                        res["R2"].append([label, mate.reference_start, mate.reference_end])
                        res["label"].add(label)
                except ValueError:
                    continue
    return res


def main(input_file: str, gtf: str, bam:str, ats_pos = 1400, mu_frag = 350, sd_frag = 50):
    u"""
    获取ATS

    按照ATS拉取范围内的transcripts。按照transcripts的范围拉取属于自己的exons

    按照ATS的范围拉取reads

    修格式导入就完了
    """

    gtf = index_gtf(gtf)

    ats = load_ats(input_file)
    print(len(ats))
    # BED("1", 11969, 14309, "+")
    for region in ats:
        iso_tbl = read_transcripts(gtf, region)

        reads = load_reads(bam, region, iso_tbl)

        if len(iso_tbl) > 0 and len(reads["label"]) > 0:
            r1_tbl = reads["R1"]
            r2_tbl = reads["R2"]

            # break all isoforms into atomic non-overlapping windows
            # exons and the transcript id
            all_wins = WinList([Window(row.start, row.end) for row in chain.from_iterable(iso_tbl.values())])
            all_win_on_gene = all_wins.split()

            # map all isoforms to the gene (atomic non-overlapping windows)
            iso_wins_list = []
            for i in iso_tbl.values():
                tmp_wins = WinList([Window(row.start, row.end) for row in i])
                tmp_wins.sort()
                iso_wins_list.append(map_iso_to_gene(tmp_wins, all_win_on_gene))

            # 需要提前将reads assign到某条转录本上
            for key in iso_tbl.keys():
                r1_wins_list = [WinList([Window(e[1], e[2]),]) for e in r1_tbl if e[0] == key]
                r2_wins_list = [WinList([Window(e[1], e[2]), ]) for e in r2_tbl if e[0] == key]
                read_labels = [1 for e in r1_tbl if e[0]==1]

                res = assign_isoform(ats_pos, iso_wins_list, r1_wins_list, r2_wins_list, read_labels, mu_frag, sd_frag)
                print(f'map_rate={res[0]}')
                print(f'iso_weights={res[1]}')

            break
    

if __name__ == '__main__':
    
    main(
        "/mnt/raid64/ATS/Personal/zhangyiming/overall_ATS_cage_R.txt",
        "/mnt/raid64/Covid19_Gravida/cellranger/Homo_sapiens/genes/genes.gtf",
        "/mnt/raid64/Covid19_Gravida/apamix/bam/CG.bam"
    )
