#!/usr/bin/env python
import sys
from Bio import SeqIO
import fftscan
import numpy as np
import pandas as pd
from time import perf_counter
import matplotlib.pyplot as plt

def allprint(f):
    print(f)
    return f

if __name__ == "__main__":
    from itertools import combinations

    if len(sys.argv) < 3:
        print("Usage: python script.py <file1.fasta> <file2.fasta> ...")
        sys.exit(1)

    allseq = [
        (record.id, record.description, record.seq)
        for f in sys.argv[1:]
        for record in SeqIO.parse(allprint(f), "fasta")
    ]

    results = []
    for a, b in combinations(allseq, 2):
        label = f"{a[0]} v {b[0]}"
        start_time = perf_counter()
        x, y = fftscan.run(a[2], b[2])
        elapsed = perf_counter() - start_time
        results.append((label, (x, y), elapsed))

    for id, desc, seq in allseq:
        print(f"{id} ({len(seq)}): {desc}")
    print("Done\n")

    labels = []
    closeness_scores = []
    lengths = []
    min_vals = []
    var_vals = []
    times = []

    for label, (x, y), elapsed in results:
        y = np.array(y)
        n_points = len(y)
        if n_points == 0:
            print(f"{label}: No data points.")
            continue
        factor = (n_points // 2) + 1
        max_val = np.max(y)
        min_val = np.min(y)
        var_val = np.var(y)
        print(f"{label} (length: {factor}):\n  Max/halfN = {max_val}\n  Min/halfN = {min_val}\n  Var/halfN = {var_val}")
        print(f"  Execution time: {elapsed:.6f} seconds\n")

        labels.append(label)
        lengths.append(factor)
        min_vals.append(min_val)
        var_vals.append(var_val)
        times.append(elapsed)
        closeness_scores.append(max_val)

    df = pd.DataFrame({
        'Label': labels,
        'Length': lengths,
        'Min': min_vals,
        'Variance': var_vals,
        'Closeness Score': closeness_scores,
        'Execution Time': times,
    })

    df.to_csv('results_summary.csv', index=False)


    sorted_data = sorted(zip(labels, closeness_scores), key=lambda x: x[1], reverse=False)

    plt.figure(figsize=(10, 6))
    plt.barh([lbl for lbl, _ in sorted_data], [score for _, score in sorted_data], color='skyblue')
    plt.xlabel("Closeness Score (Max of y)")
    plt.title("Pairwise Sequence Closeness (Higher = Closer)")
    plt.tight_layout()
    plt.show()
