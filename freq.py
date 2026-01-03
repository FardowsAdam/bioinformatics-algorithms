'''
Name : Fradows Adam 
Id : 4313204
Email : 4313204@upm.edu.sa
OS: Windows
'''


# -----------------------------------------------------------
# Program: Frequent k-mer Finder
# Description:
# This program reads DNA sequence data from a file. The file
# contains a DNA sequence (multiple lines) and a number k on 
# the last line. The program finds and prints the most frequent 
# k-length substrings (k-mers) in the sequence.
# -----------------------------------------------------------



def read_content_from_file(file_path):
    """
    Reads DNA sequence text and integer k from a given file.
    Assumes all lines except the last are text, and the last line is k.

    Args:
        file_path (str): Path to the input file.

    Returns:
        tuple: (text (str), k (int))
    """
    with open(file_path, 'r') as file:       # Open file in read mode
        content = file.read()                # Read the entire file content as a string

    content_as_list = content.strip().split('\n')  # Split file content into lines
    text = ''.join(content_as_list[:-1])           # Join all lines except the last into one string (DNA sequence)
    k = int(content_as_list[-1])                   # Convert the last line into integer k
    return text, k                                 # Return DNA sequence and k value


def most_frequent_kmers(text, k):
    """
    Finds the most frequent k-mers in the given text.

    Args:
        text (str): Input DNA sequence.
        k (int): Length of the k-mer.

    Returns:
        list: List of most frequent k-mers.
    """
    frequancy_dict = {}                # Dictionary to store k-mer counts

    # Slide a window of size k over the text
    for i in range(len(text) - k + 1):
        pattarn = text[i:i+k]          # Extract k-mer (substring of length k)
        if pattarn not in frequancy_dict:
            frequancy_dict[pattarn] = 1    # If new, set count to 1
        else:
            frequancy_dict[pattarn] += 1  # Otherwise increment count

    max_count = max(list(frequancy_dict.values()))  # Find the highest count

    # Collect all k-mers that have the highest frequency
    most_frequant =[]
    for item in frequancy_dict:
        if frequancy_dict[item] == max_count:
            most_frequant.append(item)

    return most_frequant


# ---------- Running The program ----------#

# Dataset 1: find most frequent k-mers
text1, k1 = read_content_from_file("dataset1.txt")       # Read dataset1
most_frequant1 = most_frequent_kmers(text1, k1)          # Get most frequent k-mers
print(f'The most frequant k-mer in the first file: \n {most_frequant1}')

print()  # Print blank line for separation

# Dataset 2: find most frequent k-mers
text2, k2 = read_content_from_file("dataset2.txt")       # Read dataset2
most_frequant2 = most_frequent_kmers(text2, k2)          # Get most frequent k-mers
print(f'The most frequant k-mer in the second file: \n {most_frequant2}')


most_frequant = most_frequent_kmers('CGCCTAAATAGCCTCGCGGAGCCTTATGTCATACTCGTCCT', 3)

print(most_frequant)


