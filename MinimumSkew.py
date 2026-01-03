#############################################################
# Part 2a: Finding OriC using Minimum Skew
#############################################################

# Read genome file (E. coli genome in this case)
with open("E_coli.txt", 'r') as file:
    genome = file.read().strip()


def find_min_skew_pos(genome):
    """
    Find positions in the genome where the skew diagram attains its minimum.
    Skew(i) = (# of G's - # of C's) in genome[0..i]

    Args:
        genome (str): DNA string.

    Returns:
        list: All positions where skew is minimized (OriC candidates).
    """
    n = len(genome)
    counts = [0] * (n+1)   # Initialize skew counts (0 at position 0)

    # Build skew array
    for i in range(1, n+1):
        if genome[i-1] == 'G':
            counts[i] = counts[i-1] + 1
        elif genome[i-1] == 'C':
            counts[i] = counts[i-1] - 1
        else:
            counts[i] = counts[i-1]

    min_count = min(counts)   # Find minimum skew value

    # Collect all positions where skew equals minimum
    Oric_positions = []
    for i in range(len(counts)):
        if counts[i] == min_count:
            Oric_positions.append(i)

    return Oric_positions


# Print candidate OriC positions
print("OriC candidate positions:", find_min_skew_pos(genome))



#############################################################
# Part 2b: Frequent k-mers in OriC region
#############################################################

def most_frequent_kmers(text, k):
    """
    Find the most frequent k-mers in the given text (exact match only).

    Args:
        text (str): Input DNA sequence.
        k (int): Length of k-mer.

    Returns:
        list: List of most frequent k-mers.
    """
    frequancy_dict = {}

    # Count exact k-mer occurrences
    for i in range(len(text) - k + 1):
        pattarn = text[i:i+k]
        if pattarn not in frequancy_dict:
            frequancy_dict[pattarn] = 1
        else:
            frequancy_dict[pattarn] += 1

    max_count = max(frequancy_dict.values())

    # Collect all most frequent k-mers
    most_frequant = []
    for item in frequancy_dict:
        if frequancy_dict[item] == max_count:
            most_frequant.append(item)

    return most_frequant


# Use OriC positions and extract region (500 bases after OriC candidate)
position = find_min_skew_pos(genome)
ori = genome[position[0]:position[0] + 500]

k = 9
most_freq = most_frequent_kmers(ori, k)
print("Most frequent k-mers in OriC region:", most_freq)