
def kmer2str(val, k):
    """ Transform a kmer integer into a its string representation
    :param int val: An integer representation of a kmer
    :param int k: The number of nucleotides involved into the kmer.
    :return str: The kmer string formatted
    """
    letters = ['A', 'C', 'T', 'G']
    str_val = []
    for i in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2

    str_val.reverse()
    return "".join(str_val)

# A = 00, C = 01, T = 10, G = 11
let_val = {"A": 0, "C": 1, "T": 2, "G": 3}
rev_let_val = {"A": 2, "C": 3, "T": 0, "G": 1}

def stream_kmers(text, k):
    kmer = 0
    rkmer = 0   
    for letter in text[:k-1]:
        # Forward kmer
        kmer <<= 2
        kmer += let_val[letter]
        rkmer >>= 2
        rev_letter_value = rev_let_val[letter]
        rkmer += rev_letter_value << (2 * (k - 1))
        
    mask = (1 << (2 * k)) - 1    
    for letter in text[k-1:]:
        kmer <<= 2
        kmer += let_val[letter]
        kmer &= mask
        # Reverse kmer
        rkmer >>= 2
        rev_letter_value = rev_let_val[letter]
        rkmer += rev_letter_value << (2 * (k - 1))

        yield kmer, rkmer
