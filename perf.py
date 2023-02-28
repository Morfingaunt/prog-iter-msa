#!/usr/bin/python3.10
import time

from skbio import Sequence
from skbio.sequence.distance import kmer_distance

from msa import (build_tree, get_alignment, iter_alg, kmer_distance,
                 load_datasets, prog_alg)

INDICATOR = 0.3


def test_data(msa, test_seqs):
    test_list = test_seqs
    hits = 0
    for test in test_list:
        count_for_one_seq = 0
        for seq in msa:
            al1, al2 = get_alignment(test, seq)
            if kmer_distance(Sequence(al1), Sequence(al2)) < INDICATOR:
                count_for_one_seq += 1
        pers_for_one_seq = count_for_one_seq / len(msa) * 100
        if pers_for_one_seq > 0.5:
            hits += 1
    print("Similarity: " + str("%.2f" % (hits / len(test_list) * 100)) + "%")
    print(f"Recognized {hits} out of {len(test_list)}.")


def perf_iter(filename):
    seqs = load_datasets(filename)
    tree = build_tree(seqs)
    t_start = time.perf_counter()
    msa = iter_alg(seqs, tree)
    all_time = time.perf_counter() - t_start
    print(f"iter msa: {all_time}")
    return msa


def perf_prog(filename):
    seqs = load_datasets(filename)
    tree = build_tree(seqs)
    t_start = time.perf_counter()
    msa = prog_alg(seqs, tree)
    all_time = time.perf_counter() - t_start
    print(f"prog msa: {all_time}")
    return msa


if __name__ == "__main__":
    # seqs = load_datasets("sqllong.txt")
    # build_tree(seqs, True)
    # exit()
    msa_prog = perf_prog("learn.txt")
    msa_iter = perf_iter("learn.txt")

    seqs = load_datasets("libinjection-bypasses.txt")
    print("\n======libinjection-bypasses======")
    print("   Iterative")
    test_data(msa_iter, seqs)
    print("   Progressive")
    test_data(msa_prog, seqs)

    seqs = load_datasets("sqlshort.txt")
    print("\n=============sqlshort============")
    print("   Iterative")
    test_data(msa_iter, seqs)
    print("   Progressive")
    test_data(msa_prog, seqs)

    seqs = load_datasets("sqllong.txt")
    print("\n=============sqllong=============")
    print("   Iterative")
    test_data(msa_iter, seqs)
    print("   Progressive")
    test_data(msa_prog, seqs)

    seqs = load_datasets("fake.txt")
    print("\n===============fake==============")
    print("   Iterative")
    test_data(msa_iter, seqs)
    print("   Progressive")
    test_data(msa_prog, seqs)
