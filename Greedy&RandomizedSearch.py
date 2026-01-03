"""
Name : Fradows Adam
ID   : 4313204
Email: 4313204@upm.edu.sa
OS   : Windows
Description:
This program finds DNA motifs using:
1. Greedy Motif Search
2. Randomized Motif Search (Monte Carlo)
"""

import random
import math

# ---------------------- Helper Functions ----------------------

def hamming_distance(s1, s2):
    """Count how many positions are different."""
    count = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            count += 1
    return count


def consensus(motifs):
    """Find the consensus string from motifs."""
    k = len(motifs[0])
    result = ""
    for j in range(k):
        col = [motif[j] for motif in motifs]
        counts = {'A': col.count('A'), 'C': col.count('C'),
                  'G': col.count('G'), 'T': col.count('T')}
        result += max(counts, key=counts.get)
    return result


def score(motifs):
    """Return total distance of motifs from consensus."""
    cons = consensus(motifs)
    total = 0
    for motif in motifs:
        total += hamming_distance(cons, motif)
    return total


def build_profile(motifs):
    """Build profile matrix with pseudocounts."""
    k = len(motifs[0])
    t = len(motifs)
    profile = [[1] * k for _ in range(4)]  # pseudocounts
    letters = ['A', 'C', 'G', 'T']

    for j in range(k):
        for motif in motifs:
            ch = motif[j]
            profile[letters.index(ch)][j] += 1

    for i in range(4):
        for j in range(k):
            profile[i][j] /= (t + 4)
    return profile


# ---------------------- Core Algorithms ----------------------

def most_probable_kmer(text, k, profile):
    """Find the most probable k-mer based on profile."""
    max_prob = -1
    best = text[:k]
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        prob = 1
        for j in range(k):
            if kmer[j] == 'A':
                prob *= profile[0][j]
            elif kmer[j] == 'C':
                prob *= profile[1][j]
            elif kmer[j] == 'G':
                prob *= profile[2][j]
            elif kmer[j] == 'T':
                prob *= profile[3][j]
        if prob > max_prob:
            max_prob = prob
            best = kmer
    return best


def greedy_motif_search(dna_list, k, t):
    """Run Greedy Motif Search."""
    best_motifs = [dna[:k] for dna in dna_list]
    n = len(dna_list[0])

    for i in range(n - k + 1):
        motifs = [dna_list[0][i:i+k]]
        for j in range(1, t):
            profile = build_profile(motifs)
            motifs.append(most_probable_kmer(dna_list[j], k, profile))
        if score(motifs) <= score(best_motifs):
            best_motifs = motifs[:]
    return best_motifs


def randomized_motif_search(dna_list, k, t):
    """Run Randomized Motif Search once."""
    n = len(dna_list[0])
    positions = random.choices(range(n - k + 1), k=t)
    motifs = [dna_list[i][positions[i]:positions[i]+k] for i in range(t)]
    best = motifs[:]
    best_score = score(best)

    while True:
        profile = build_profile(motifs)
        motifs = [most_probable_kmer(dna_list[i], k, profile) for i in range(t)]
        if score(motifs) < best_score:
            best = motifs[:]
            best_score = score(motifs)
        else:
            return best, best_score


def randomized_search_monte_carlo(dna_list, k, t, runs=200):
    """Run Randomized Motif Search many times."""
    best = None
    best_score = math.inf
    for _ in range(runs):
        motifs, s = randomized_motif_search(dna_list, k, t)
        if s < best_score:
            best = motifs
            best_score = s
    return best, best_score


# ---------------------- Example Run ----------------------

if __name__ == "__main__":
    dna = "GGCGTTCAGGCA AAGAATCAGTCA CAAGGAGTTCGC CACGTCAATCAC CAATAATATTCG"
    dna_list = dna.split()
    k = 3
    t = 5

    print("\n--- Greedy Motif Search ---")
    best_greedy = greedy_motif_search(dna_list, k, t)
    print("Best Motifs (Greedy):", best_greedy)

    print("\n--- Randomized Motif Search ---")
    best_random, score_val = randomized_search_monte_carlo(dna_list, k, t)
    print("Best Motifs (Randomized):", best_random)
