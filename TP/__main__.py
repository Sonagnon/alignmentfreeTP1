from TP.loading import load_directory, load_fasta
from TP.kmers import stream_kmers, kmer2str, compare_kmers, encode_kmer, encode_nucl, enumerate_kmers, compare_kmers, create_index, compare_sequences




def jaccard(fileA, fileB, k):
    '''Fonction qui calcule la distance de Jaccard
    Entrée : elle prend en entrée les deux séquence à comparer et la taille du kmers
    Sortie: elle retourne la valeur j de la distance calculée'''
    j = 0

    intersect_value, union_value = compare_sequences(fileA, fileB, k)
    
    j = (intersect_value / union_value)
        
    return j




if __name__ == "__main__":
    print("Computation of Jaccard similarity between files")

    # Load all the files in a dictionary
    files = load_directory("data")
    k = 21
    
    print("Computing Jaccard similarity for all pairs of samples")
    filenames = list(files.keys())
    #print(files)
    print(filenames)
    sequence1 = []
    sequence2 = []
    for i in range(len(files)):
        for j in range(i+1, len(files)):
            sequence1.append(files[filenames[i]])
            sequence2.append(files[filenames[j]])
            
            #j_dist = jaccard(files[filenames[i]], files[filenames[j]], k)
            #j_dist = jaccard(sequence1, sequence2, k)
            n, b = compare_sequences(sequence1, sequence2, k)
            print(filenames[i], filenames[j], n, b)
            #print(filenames[i], filenames[j])
