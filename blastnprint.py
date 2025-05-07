#!/usr/bin/env python
import os
import sys
import subprocess
import tempfile
from itertools import combinations
from time import perf_counter

from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt

def run_blastn(query_seq, subject_seq):
    with tempfile.TemporaryDirectory() as tmpdir:
        query_file = os.path.join(tmpdir, "query.fasta")
        subject_file = os.path.join(tmpdir, "subject.fasta")
        db_prefix = os.path.join(tmpdir, "subject_db")
        output_file = os.path.join(tmpdir, "blast_result.txt")

        with open(query_file, "w") as qf:
            qf.write(f">query\n{query_seq}\n")
        with open(subject_file, "w") as sf:
            sf.write(f">subject\n{subject_seq}\n")
        subprocess.run(["makeblastdb", "-in", subject_file, "-dbtype", "nucl", "-out", db_prefix],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            "blastn", "-query", query_file, "-db", db_prefix,
            "-outfmt", "6 bitscore", "-out", output_file
        ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        try:
            with open(output_file) as f:
                scores = [float(line.strip()) for line in f if line.strip()]
            return max(scores) if scores else 0.0
        except FileNotFoundError:
            return 0.0

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <file1.fasta> <file2.fasta> ...")
        sys.exit(1)

    allseq = [
        (record.id, record.description, str(record.seq))
        for f in sys.argv[1:]
        for record in SeqIO.parse(f, "fasta")
    ]

    results = []
    for a, b in combinations(allseq, 2):
        label = f"{a[0]} v {b[0]}"
        start_time = perf_counter()
        score = run_blastn(a[2], b[2])
        elapsed = perf_counter() - start_time
        results.append((label, score, elapsed))

    for id, desc, seq in allseq:
        print(f"{id} ({len(seq)}): {desc}")
    print("Done\n")

    df = pd.DataFrame(results, columns=["Label", "Score", "Time"])

    df.to_csv('results_summary_blast.csv', index=False)

    sorted_data = df.sort_values(by="Score", ascending=False)

    plt.figure(figsize=(10, 6))
    plt.barh(sorted_data["Label"], sorted_data["Score"], color='skyblue')
    plt.xlabel("Closeness Score (Max Bit Score from BLASTN)")
    plt.title("Pairwise Sequence Closeness via BLASTN")
    plt.tight_layout()
    plt.show()