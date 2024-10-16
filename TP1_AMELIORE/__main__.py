from TP1_AMELIORE.loading import load_directory, load_fasta
from TP1_AMELIORE.kmers import compare_sequences, enumerate_kmers, create_index, encode_kmer_binary, generate_canonical_kmers, canonical_kmer, reverse_complement, encode_nucl

def calculate_jaccard(sequence1, sequence2, k):
    '''
    Fonction qui calcule la distance de jaccard
    '''

    inter, uni = compare_sequences(sequence1, sequence2, k)

    return inter / uni if uni != 0 else 0






if __name__ == "__main__":
    print("Computation of Jaccard similarity between files")
    
    # Load all the files in a dictionary
    files = load_directory("data")
    k = 21
    filenames = list(files.keys()) 

    print("Computing Jaccard similarity for all pairs of samples")
    for i in range(len(filenames)):
        sequence1 = files[filenames[i]]
        for j in range(i, len(filenames)):
            sequence2 = files[filenames[j]]
            jaccard_dist = calculate_jaccard(sequence1, sequence2, k)
            print(filenames[i], filenames[j], jaccard_dist)
