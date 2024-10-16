def encode_nucl(letter):
    '''Fonction qui encode les nucléotides en binaire  '''
    encoding = {'A': 0b00,  'C': 0b01,  'T': 0b10,  'G': 0b11  }
    return encoding[letter]

def reverse_complement(kmer):
    ''' Fonction qui génère le complément inversé d'un kmer '''
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(kmer))

def canonical_kmer(kmer):
    '''Fonction qui va nous sortir le k-mer canonique entre le kmer et son reserve complement
    sachant que le kmer canonique sera le kmer qui sera le plus petit lexicographiquement parlant.'''
    
    rev_comp = reverse_complement(kmer)
    return min(kmer, rev_comp)  
        

def generate_canonical_kmers(sequence, k):
    ''' Fonction qui va nous générer tous les k-mers canoniques d'une séquence '''
    #kmers = []
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        canonical = canonical_kmer(kmer)
        #kmers.append(canonical)
        yield canonical



def encode_kmer_binary(kmer):
    ''' Fonction qui encode un k-mer en binaire en utilisant la fonction encode_nucl()'''
    binary_kmer = 0 
    binary_string = ''
    for letter in kmer:
        binary_kmer = (binary_kmer << 2) | encode_nucl(letter)  # Décalage de 2 bits à chaque étape
        binary_string += format(binary_kmer & 0b11, '02b')  # '02b' signifie 2 bits, avec des zéros à gauche si nécessaire

    return binary_string 



def create_index(seq, k):
    """ Fonction qui Cré un index des k-mers encodés d'une séquence donnée.
    
    :param seq: La séquence d'entrée.
    :param k: La longueur des k-mers.
    :return: Un dictionnaire avec les k-mers encodés comme clés et leurs fréquences comme valeurs.
    Nous avons fait le choix d'encoder aussi les kmers dans l'index, pour faciliter la comparaison des kmers dans la fonction suivantes.
    """
    index = {}
    #kmers = generate_canonical_kmers(seq, k)
    # Itérer sur tous les k-mers dans la séquence
    for kmer in generate_canonical_kmers(seq, k):
        #kmers = generate_canonical_kmers(seq[i:i + k], k)  # ici on génère le k-mer canonique
        the_kmer = encode_kmer_binary(kmer) #j'essaie d'encoder le kmer avant de le compter
        # je met à jour la fréquence du k-mer dans l'index 
        if the_kmer in index:
            index[the_kmer] += 1
        else:
            index[the_kmer] = 1
            
    return index

def enumerate_kmers(seq, k):
    '''Fonction qui retourne tous les k-mers d'une séquence ou un text 
        Entrée : la séquence et la taille k des kmers
        Sortie: un générateur contenant les k-mers
    '''
    #mask = (1 << (2 * k)) - 1 
    for kmer in generate_canonical_kmers(seq, k):
        the_kmer = kmer
        kmer_binaire = encode_kmer_binary(the_kmer)
        yield kmer_binaire


def compare_sequences(seq1, seq2, k):
    ''' Fonction qui compare deux séquences fournie en se basant sur la méthodes des k-mers.
    Entrées: seq1 et seq2 les deux séquences à comparer
                k, la taille des kmers
    Sortie : on renvoit
            intersect , la taille de l'intersection entre les deuc séquence
            union, la taille de l'union des deux séquences
    '''
    index = create_index(seq1, k)  
    
    intersect = 0
    union = sum(index.values()) 

    for kmer in enumerate_kmers(seq2, k): 
        if kmer in index and index[kmer] > 0:
            intersect += 1
            index[kmer] -= 1
        else:
            union += 1

    return intersect,  union
