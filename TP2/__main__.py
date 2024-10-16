from TP2.loading import load_directory, load_fasta
from TP2.kmers import compare_sketches, create_sketch, encode_kmer_binary, generate_canonical_kmers, canonical_kmer, reverse_complement, encode_nucl, hash_binary_kmer


def calculate_jaccard_2(sequence1, sequence2, k, sketch_size):
    '''Fonction qui Calcule la distance de Jaccard entre les sketches de deux séquences'''
    sketch1 = create_sketch(sequence1, k, sketch_size)
    sketch2 = create_sketch(sequence2, k, sketch_size)
    inter, uni = compare_sketches(sketch1, sketch2)
    return inter / uni if uni != 0 else 0.0





if __name__ == "__main__":
    print("Computation of Jaccard similarity between files")
    
    # Load all the files in a dictionary
    files = load_directory("data")
    k = 21
    sketch_size = 500
    filenames = list(files.keys()) 
    print(filenames)
    print("Computing Jaccard similarity for all pairs of samples")
    filenames = list(files.keys())
    for i in range(len(files)):
        for j in range(i, len(files)):
            jaccard_dist = calculate_jaccard_2(files[filenames[i]], files[filenames[j]], k, sketch_size)
            print(filenames[i], filenames[j], jaccard_dist)
            
            

