
def kmer2str(val, k):
    """ Transform a kmer integer into a its string representation
    :param int val: An integer representation of a kmer
    :param int k: The number of nucleotides involved into the kmer.
    :return str: The kmer string formatted
    """
    letters = ['A', 'C', 'T', 'G']
    str_val = []
    for _ in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2

    str_val.reverse()
    return "".join(str_val)


def stream_kmers(text, k):
    # --- To complete ---
    pass
    
def compare_kmers(km1, km2):
    for i in range(len(km1)):
        if km1[i] != km2[i]:
            return False
    return True


def encode_kmer(kmer):
    '''Fnction qui encode les  k-mers d'une séquence données
        Entrée: la sequence seq et la taille k des k-mers
        Sortie: un objet de type generator
    '''
    kmer_enc = 0
    #for letter in seq[:k]:  
    for letter in kmer:
        kmer_enc <<= 2
        kmer_enc += encode_nucl(letter)
    return kmer_enc


def encode_nucl(letter):
    encoding = {'A': 0b00,  'C': 0b01,  'T': 0b10,  'G': 0b11  }
    return encoding[letter]

def enumerate_kmers(seq, k):
    '''Fonction qui retourne tous les k-mers d'une séquence ou un text
        Entrée : la séquence et la taille k des kmers
        Sortie: un générateur contenant les k-mers
    '''
    mask = (1 << (2 * k)) - 1
    #kmer = encode_kmer(seq[0:k-1], k)
    kmer = encode_kmer(seq, k)
    yield kmer  # Premier k-mer
    for i in range(1, len(seq) - k + 1):
        kmer = (kmer << 2) & mask  # Décaler à gauche et appliquer le masque
        kmer += encode_nucl(seq[i + k - 1])  # Ajouter le nouveau nucléotide
        yield kmer


def compare_kmers(km1, km2):
    '''Fonction qui compare deux k-mers entre eux
        km1 : kmer 1
        km2 : kmer 2

        Sortie: Vrai ou Faux
    '''
    for i in range(len(km1)):
        if km1[i] != km2[i]:
            return False
    return True




def create_index(seq, k):
    """ Fonction qui Cré un index des k-mers encodés d'une séquence donnée.
    :param seq: La séquence d'entrée.
    :param k: La longueur des k-mers.
    :return: Un dictionnaire avec les k-mers encodés comme clés et leurs fréquences comme valeurs.

    Nous avons fait le choix d'encoder aussi les kmers dans l'index, pour faciliter la comparaison des kmers dans la fonction suivantes.
    """
    index = {}
    # Itérer sur tous les k-mers dans la séquence
    for i in range(len(seq) - k + 1):
        kmer = encode_kmer(seq[i:i + k], k)  # Encoder le k-mer

        # Ajouter ou mettre à jour la fréquence du k-mer dans l'index
        if kmer in index:
            index[kmer] += 1
        else:
            index[kmer] = 1
    return index



def compare_sequences(seq1, seq2, k):
    ''' Fonction qui compare deux séquences fournie en se basant sur la méthodes des k-mers.
    Entrées: seq1 et seq2 les deux séquences à comparer
                k, la taille des kmers
    Sortie : on renvoit
            intersect , la taille de l'intersection entre les deuc séquence
            union, la taille de l'union des deux séquences
    '''
    index = create_index(seq1, k)
    print(index)
    intersect = 0
    union = sum(index.values())
    for kmer in enumerate_kmers(seq2, k):

        if kmer in index and index[kmer] > 0:
            intersect += 1
            index[kmer] -= 1
        else:
            union += 1

    return intersect,  union
