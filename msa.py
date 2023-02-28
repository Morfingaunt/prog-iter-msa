#!/usr/bin/python3.10
import time
from functools import partial

from Bio.Align import PairwiseAligner
from iab.algorithms import progressive_msa, progressive_msa_and_tree
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import average, dendrogram
from skbio import DNA, DistanceMatrix, Sequence, TabularMSA, TreeNode
from skbio.sequence.distance import kmer_distance

# k - blocksize
kmer_distance = partial(kmer_distance, k=3, overlap=False)


def erase_dash(aln1, aln2, max_len):
    align1 = align2 = ""
    pairs = list(zip(aln1, aln2))
    while len(pairs) > max_len:
        if ("-", "-") in pairs:
            pairs.remove(("-", "-"))
        else:
            pairs.pop()
    align1 = ''.join(list(map(lambda tup: tup[0], pairs)))
    align2 = ''.join(list(map(lambda tup: tup[1], pairs)))
    return align1[:max_len], align2[:max_len]


def get_alignment(aln1, aln2, max_len=0):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    alignnment = aligner.align(str(aln1), str(aln2))[0]
    alignnments = str(alignnment).split("\n")
    if max_len > 0:
        return erase_dash(alignnments[0], alignnments[2], max_len)
    return alignnments[0], alignnments[2]


def nw_align(aln1: TabularMSA | Sequence, aln2: TabularMSA | Sequence):
    max_len1 = 0
    max_len2 = 0

    if isinstance(aln1, DNA | Sequence):
        aln1 = TabularMSA([aln1], metadata=aln1.metadata)
    elif isinstance(aln1, TabularMSA):
        max_len1 = len(aln1[0])
    if isinstance(aln2, DNA | Sequence):
        aln2 = TabularMSA([aln2], metadata=aln2.metadata)
    elif isinstance(aln2, TabularMSA):
        max_len2 = len(aln2[0])
    max_len = max(max_len1, max_len2)
    msa_left = []
    msa_right = []
    msa_left.extend(aln1)
    msa_right.extend(aln2)

    for idx, val_right in enumerate(msa_right):
        for jdx, val_left in enumerate(msa_left):
            al1, al2 = get_alignment(val_left, val_right, max_len)
            msa_left[jdx] = \
                aln1[jdx]._constructor(
                    sequence=al1,
                    metadata=aln1[jdx].metadata,
                    positional_metadata=None
                )
            msa_right[idx] = \
                aln2[idx]._constructor(
                    sequence=al2,
                    metadata=aln2[idx].metadata,
                    positional_metadata=None
                )
    msa_left.extend(msa_right)
    return TabularMSA(msa_left), None, None


def iter_alg(seq_list, guide_tree=None, iterations=3):
    tree = guide_tree
    for i in range(iterations):
        msa, tree = progressive_msa_and_tree(
            seq_list,
            pairwise_aligner=nw_align,
            metric=kmer_distance,
            guide_tree=tree
        )
    return msa


def prog_alg(seq_list, guide_tree=None, save_dm=False):
    msa = progressive_msa(
        seq_list,
        pairwise_aligner=nw_align,
        metric=kmer_distance,
        guide_tree=guide_tree
    )
    if save_dm:
        msa_dm = DistanceMatrix.from_iterable(
            msa,
            metric=kmer_distance,
            key="id"
        )
        msa_dm.plot(cmap="Greens")
        plt.savefig('fooprog.png')
    return msa


def build_tree(seqs, save_tree=False):
    # t_start = time.perf_counter()
    guide_dm = DistanceMatrix.from_iterable(
        seqs,
        metric=kmer_distance,
        key='id'
    )
    guide_lm = average(guide_dm.condensed_form())
    guide_tree = TreeNode.from_linkage_matrix(guide_lm, guide_dm.ids)
    if save_tree:
        dendrogram(guide_lm, labels=guide_dm.ids, orientation='right',
                      link_color_func=lambda x: 'black')
        plt.savefig('tree.png')
    # print(guide_tree.ascii_art())
    # all_time = time.perf_counter() - t_start
    # print(f"tree: {all_time}")
    return guide_tree


def load_datasets(filename):
    with open(filename, "r") as file:
        data_sets = file.readlines()
        seqs = []
        for idx, value in enumerate(data_sets):
            value = value.rstrip('\n')
            seqs.append(Sequence(value, {"id": f"s{idx}"}))
        return seqs
