from loading import load_directory
from kmers import stream_kmers
from collections import Counter

def jaccard(fileA, fileB, k):
    j = 0
    
    # Store k-mers of fileA in a dictionary
    index_fileA = index_file(fileA, k)

    # Find k-mers unique for fileA, unique for fileB and in intersection of the two
    A_only, A_and_B, B_only = intersect_index(index_fileA, fileB, k)

    # Jaccard index: intersection/ (|A| + |B| - intersection)
    j = A_and_B / (A_only + A_and_B + B_only)

    return j

def list_kmers(sequences, k):
    kmers = []
    for seq in sequences:
        kmers.extend([min(kmer, rkmer) for kmer, rkmer in stream_kmers(seq, k)])
    return kmers

def index_file(sequences, k):
    return Counter(list_kmers(sequences, k))

def intersect_index(index, sequences, k):
    """ Create an index containing all the kmers of the input sequences
    """
    index_uniq = index.copy()
    query_uniq = 0
    intersection = 0

    for seq in sequences:
        for kmer, rkmer in stream_kmers(seq, k):
            minmer = min(kmer, rkmer)

            # Query not in index
            if minmer not in index_uniq:
                query_uniq += 1
            # Query in index => intersection
            else:
                intersection += 1
                index_uniq[minmer] -= 1
                if index_uniq[minmer] == 0:
                    del index_uniq[minmer]

    return sum(index_uniq.values()), intersection, query_uniq

if __name__ == "__main__":
    # Load all the files in a dictionary
    files = load_directory("data")
    k = 21
    
    filenames = list(files.keys())
    for i in range(len(files)):
        for j in range(i+1, len(files)):       
            jac = jaccard(files[filenames[i]], files[filenames[j]], k)
            print(filenames[i], filenames[j], jac)
