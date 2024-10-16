''''Je vais réutiliser mes anciens codes qui me servaient à créer les kmers puis je 
 redéfinir/réadapter certaines anciennes focntions'''''

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

def hash_binary_kmer(binary_kmer):
    '''Fonction qui va hacher mon k-mer sur 16bit qui est déjà sous forme binaire en utilisant SHA-1 pour une distribution uniforme'''
    return int(hashlib.sha1(binary_kmer.encode('utf-8')).hexdigest(), 16)



def create_sketch(sequence, k, sketch_size):
    '''Cette fonction va creer un sketch en filtrant les kmers de la séquence et en conservant les plus petite valeurs de kmers haché'''
    sketch = []
    
    for kmer in generate_canonical_kmers(sequence, k):        #donc je commence par générer les kmers canoniques à partir de la séquence
        #et pour chaque kmer, je l'encode je le hache ensuite et j'utilise heapq pour l'ajouter dans ma liste sketch ou l'enlever  selon si il est petit ou pas 
        binary_kmer = encode_kmer_binary(kmer)  
        hashed_kmer = hash_binary_kmer(binary_kmer) 
        
        if len(sketch) < sketch_size:
            heapq.heappush(sketch, -hashed_kmer)  
        else:
            heapq.heappushpop(sketch, -hashed_kmer)
    
    return sorted([-x for x in sketch]) 


def compare_sketches(sketch1, sketch2):
    '''Cette fonction va comparer les deux sketches et renvoyer l'intersection et l'union '''
    
    i, j = 0, 0
    intersect = 0
    union = 0
    
    while i < len(sketch1) and j < len(sketch2):
        if sketch1[i] == sketch2[j]:
            intersect += 1
            union += 1
            i += 1
            j += 1
        elif sketch1[i] < sketch2[j]:
            union += 1
            i += 1
        else:
            union += 1
            j += 1
    
    # Si il reste des éléments à la sortie de notre boucle on les ajoute dans union
    union += len(sketch1) - i + len(sketch2) - j
    
    return intersect , union 
