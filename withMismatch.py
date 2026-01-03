'''
Name : Fradows Adam 
Id : 4313204
Email : 4313204@upm.edu.sa
OS: Windows
'''

#############################################################
# Part 1: Frequent Words with Mismatches Problem
#############################################################

def HammingDistance(s1, s2, d):
    """
    Compute whether two strings differ in at most `d` positions.

    Args:
        s1 (str): First string.
        s2 (str): Second string (same length as s1).
        d (int): Maximum allowed mismatches.

    Returns:
        bool: True if Hamming distance between s1 and s2 <= d, else False.
    """
    count = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:   # mismatch found
            count += 1
    return count <= d


def most_frequent_kmers(text, k, d):
    """
    Find the most frequent k-mers in `text` allowing up to `d` mismatches.

    Args:
        text (str): Input DNA sequence.
        k (int): Length of k-mer.
        d (int): Maximum allowed mismatches.

    Returns:
        list: List of k-mers that occur most frequently (with mismatches allowed).
    """
    frequancy_dict = {}  # Dictionary to store exact k-mer counts

    # Step 1: Count exact k-mer occurrences
    for i in range(len(text) - k + 1):
        pattarn = text[i:i+k]              # Extract substring of length k
        if pattarn not in frequancy_dict:
            frequancy_dict[pattarn] = 1
        else:
            frequancy_dict[pattarn] += 1

    # Step 2: Count frequencies including mismatches
    frequancy_dict_with_mismatch = {}

    for pattarn in frequancy_dict:
        total = 0
        # Compare with every other k-mer
        for pattarn_ in frequancy_dict:
            if HammingDistance(pattarn, pattarn_, d):   # If within d mismatches
                total += frequancy_dict[pattarn_]
        frequancy_dict_with_mismatch[pattarn] = total

    # Step 3: Find the maximum frequency
    max_count = max(frequancy_dict_with_mismatch.values())

    # Step 4: Collect all k-mers with max frequency
    most_frequant = []
    for item in frequancy_dict_with_mismatch:
        if frequancy_dict_with_mismatch[item] == max_count:
            most_frequant.append(item)

    return most_frequant


# Example test run
result = most_frequent_kmers('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 1)
print("Most frequent k-mers with mismatches:", result)