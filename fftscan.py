#!/usr/bin/env python
import sys
import numpy as np
import itertools
from collections import Counter
from Bio import SeqIO
from scipy.fft import fft, ifft


#nue = "ATCG"
#encoder = {''.join(a) : i for i, a in enumerate(itertools.product(nue, repeat=4))}

nueclotide_to_signal= {
        'A': 1,
        'T': -1,
        'G': 1j,
        'C': -1j
    }



def fill_shape(A, B):
  if A.shape[0] < B.shape[0]:
      A = np.concatenate((A, np.zeros(B.shape[0] - A.shape[0])))
  elif A.shape[0] > B.shape[0]:
      B = np.concatenate((B, np.zeros(A.shape[0] - B.shape[0])))
  return A, B

def encode_nueclotide_sequence(sequenceA, sequenceB):
  A = Counter(sequenceA)
  B = Counter(sequenceB)
  T = A + B
  le = len(sequenceA) + len(sequenceB)
  # adjust for noise
  A = {k: ((0.25*v) / len(sequenceA)) * nueclotide_to_signal.get(k) for k, v in A.items()}
  B = {k: ((0.25*v) / len(sequenceB)) * nueclotide_to_signal.get(k) for k, v in B.items()}
  _A = np.fromiter(
      (A[nuclotide] for nuclotide in sequenceA), dtype=complex
  )
  _B = np.fromiter(
      (B[nuclotide] for nuclotide in sequenceB), dtype=complex
  )
  return _A, _B

def compare_fft(x, y):
    X = fft(x)
    Y = fft(y)
    return ifft(Y * np.conjugate(X)).real ,ifft(X * np.conjugate(Y)).real


# This technique requires just three FFTs, two
# forward and one inverse.

def run(seq1, seq2):

    x, y = encode_nueclotide_sequence(seq1, seq2)
    x, y = fill_shape(x, y)
    # you could also apply this method to protiens by doing it 5 times with 4 diffrent protiens each

    _y1, _y2 = compare_fft(x, y)
    _x = np.arange(-len(_y1) + 1, len(_y2))
    return _x, np.concatenate((_y1[:-1], _y2))





if __name__ == "__main__":
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.rcParams['agg.path.chunksize'] = 10000
    mpl.rcParams['path.simplify'] = True
    mpl.rcParams['path.simplify_threshold'] = 0.5
    if len(sys.argv) != 3:
        print("Usage: python script.py <file1> <file2>")
        sys.exit(1)
    file1 = sys.argv[1]
    file2 = sys.argv[2]

    for record in SeqIO.parse(file1, "fasta"):
        desc1 = record.description
        seq1 = record.seq

    for record in SeqIO.parse(file2, "fasta"):
        desc2 = record.description
        seq2 = record.seq
    print(f"{file1} size : {len(seq1)}")
    print(desc1)
    plt.title(f"{desc1} v {desc2}")
    print(f"{file2} size : {len(seq2)}")
    print(desc2)
    x,y = run(seq1, seq2)

    print("HERE")
    plt.plot(x, y)
    plt.ylim(0, max(10,max(y)))
    plt.show()









