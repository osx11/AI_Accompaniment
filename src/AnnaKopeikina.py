import random
import mido
from mido import Message ,MetaMessage
import music21

INPUT_FILE_NAME = 'input2.mid'
OUTPUT_FILE_NAME = 'output2.mid'

midi = mido.MidiFile(INPUT_FILE_NAME, clip=True)
midi.type = 1

CHORD_MINOR = [0, 3, 7]
CHORD_MAJOR = [0, 4, 7]
CHORD_DIM = [0, 3, 6]
CHORD_ALL = [CHORD_MINOR, CHORD_MAJOR, CHORD_DIM]
CHORD_DURATION = midi.ticks_per_beat * 2

score = music21.converter.parse(INPUT_FILE_NAME)
key = score.analyze('key')

if key.mode == "minor":
    MUSIC_SCALE = (key.tonic.midi+3) % 12
else:
    MUSIC_SCALE = key.tonic.midi % 12


def get_average_velocity():
    """Return the average level of sound from original song."""
    note_count=0
    avg_velocity=0
    for note in midi.tracks[1]:
        if isinstance(note, Message) and ( note.type == 'note_on'):
            avg_velocity+=note.velocity
            note_count+=1
    return int(avg_velocity/note_count)


def get_average_octave():
    """Return the average range of notes from original song."""
    note_count=0
    avg_octave=0
    for note in midi.tracks[1]:
        if isinstance(note, Message) and ( note.type == 'note_on'):
            avg_octave+=int(note.note/12)
            note_count+=1
    return int(avg_octave/note_count)


class Chord:
    def __init__(self, root_note, chord_type):
        """Class of triad chord
        
            Keyword arguments:
            root_node -- the start note for the chord
            chord_type -- arry or intends for the chords other notes
        """
        self.root_note= root_note
        self.type= chord_type
        self.note_list= [(root_note + note) % 12 for note in chord_type]
        pass
    def has_note(self, note): 
        """Return if the chord has given note in his notes """
        if(note==None):
             return False
        for chord_note in self.note_list:
            if(chord_note == note % 12):
                return True
        return False
    def __eq__(self, other): 
        """Return if this instance is equal to other chord. Overrided standard method for simpler comparasion"""
        return self.root_note == other.root_note and self.type == other.type

class Accompanement:
    """Class of accompanement, consisting of best chords for the song, size etc
        
            Keyword arguments:
            scale -- The scale for the song, needed for calculating sonsonant chords
            size -- amount of tacts, needed for calculating the final amount of chords
    """
    def __init__(self, scale, size):
        self.scale= scale
        self.size=size
        self.consonantChords = [] #The consonant chords are the appropriate chords for the music. We get them by formula, derived from the circe of fifth
        self.consonantChords.append(Chord(scale % 12, CHORD_MAJOR))
        self.consonantChords.append(Chord((scale + 5) % 12, CHORD_MAJOR))
        self.consonantChords.append(Chord((scale + 7) % 12, CHORD_MAJOR))
        self.consonantChords.append(Chord((scale + 9) % 12, CHORD_MINOR))
        self.consonantChords.append(Chord((scale + 2) % 12, CHORD_MINOR))
        self.consonantChords.append(Chord((scale + 4) % 12, CHORD_MINOR))
        self.consonantChords.append(Chord((scale + 11) % 12, CHORD_DIM))
        pass
    
    """ Return the fitting chord, which contains given note"""
    def get_consonant_chord(self,note):
        for chord in self.consonantChords:
            if(chord.has_note(note)):
                return chord
        return random.choice(self.consonantChords)

    """ Return all notes used in chords"""
    def print_all_notes(self):
        note_list = []
        for chord in self.consonantChords:
            for note in chord.note_list:
                if note not in note_list:
                    note_list.append(note)

    """ Return if the consonant chords contain the given chord"""
    def has_in_consonant_chords(self, chord):
        for existing_chord in self.consonantChords:
            if(existing_chord==chord):
                return True
        return False
        
    """ Return if the consonant chords contain the given note inside any of them"""
    def has_note_in_consonant_chords(self, note):
        if(note==None):
            return False 
        for existing_chord in self.consonantChords:
            if(existing_chord.has_note(note%12)):
                return True
        return False

""" Return the amount of half-tacks inside the music"""
def get_notes_amount(track):
    beats = 0
    for msg in track:
        if type(msg) is Message:
            beats +=msg.time
    length= (beats+CHORD_DURATION-1)//CHORD_DURATION
    return length

"""Get all the notes of the track, which are needed for computing the chords"""
def compute_beats(track):
    length=get_notes_amount(track)
    notes_list=[None]*length
    beats = 0
    last_note=0
    for msg in track:
       if type(msg) is Message:
            if(beats % CHORD_DURATION == 0 and msg.type =="note_on" ):
               if(notes_list[beats//CHORD_DURATION] is None):
                   notes_list[beats//CHORD_DURATION] = msg.note
               
            if(msg.type=="note_off"):
                last_note= msg.note
            beats +=msg.time
    if(beats % CHORD_DURATION == 0):
        notes_list[-1]=last_note
    return notes_list


"""Get the tempo of the original track"""
tempo=0
for track in midi.tracks:
    for msg in track:
        if isinstance(msg, MetaMessage) and msg.type=="set_tempo":
            tempo=msg.tempo

tracks=[MetaMessage('set_tempo',tempo=tempo,time=0)]

song_notes=compute_beats(midi.tracks[1])
velocity= int(get_average_velocity()*0.9) #The amount of how loud will be accompanement
average_displacement= 12*(get_average_octave()-1) #The range where will be notes of the accompanement

final_genom_chord = Accompanement(MUSIC_SCALE, len(song_notes))
max_note=120 


class Chromosome: 
    """Class of chromosome needed for genetic algorithm. Chromosome itself represents an accompanement, which consists of chords 
        
            Keyword arguments:
            size -- amount of genes in accompanement
    """
    def __init__(self, size):
        self.size = size
        self.genes_pool = [None]*size
        self.rating= 0
        self.generate_random_genes()

    """Generate random genes using random note in range from 0 to 120 and type of chord to get the random chord for each gene """
    def generate_random_genes(self):
        for i in range(self.size):
            rand_note  = random.randint(0,max_note)
            rand_chord = Chord(rand_note, random.choice(CHORD_ALL))
            self.genes_pool[i] = rand_chord

    """Overrated method of equals to other chromosome for quicker comparison"""
    def __eq__(self, other):
        return self.rating == other.rating 

    """Overrated method of less than to other chromosome for quicker comparison"""
    def __lt__(self, other):
        return self.rating < other.rating

"""Return the population of chromosomes with fixed amount of genes"""
def create_population(population_size, chromosome_size):
    population = [None]*population_size
    for i in range(population_size):
        population[i]= Chromosome(chromosome_size);
    return population

"""Calculate rating for each of the chromosomes. This is the fitness function of the generic algorithm"""
def calculate_rating(population):
    for chromosome in population:
        chromosome.rating= chromosome.size
        for i in range(chromosome.size):
            if final_genom_chord.has_in_consonant_chords(chromosome.genes_pool[i]): #if this gene is in consonant fitting chords, we reward it
                chromosome.rating-=0.5
                if(song_notes[i] is None):
                    chromosome.rating-=0.5
                    continue
                if(not final_genom_chord.has_note_in_consonant_chords(song_notes[i]) or chromosome.genes_pool[i].has_note(song_notes[i])): #if this genes chord is fitting to the music, we reward it
                    chromosome.rating-=0.5

"""The constrains for the fitness algorithm"""    
population_size=128
survivors = [None] * (population_size // 4) #the size of the best elite chromosomes

population= create_population(population_size, final_genom_chord.size) #Creating the original population

def select(population,survivors):
    """Select the top first chromosomes from population array""" 
    size= len(survivors)
    for i in range(size):
        survivors[i]=population[i]

def get_parent_index(parents, other_parent_index):
    """Get the indes of random chromosome from parents. Other parent index is needed to having two different parents""" 
    size= len(parents)
    while True: #We use random and while true instead of "for i in index" because we need the most uniquecombinations while crossing
        index = random.randint(0, size-1)
        if other_parent_index is None or other_parent_index!=index:
            return index

def cross(chromosome1, chromosome2):
    """Return the child from two cromosomes formed by slicing the parents and given the child the parts of two parents chromosomes. This is crossingover of algorithm""" 
    size = chromosome1.size
    point= random.randint(0,size-1)
    child = Chromosome(size)
    for i in range(point):
        child.genes_pool[i]=chromosome1.genes_pool[i]
    for i in range(point,size):
        child.genes_pool[i]=chromosome2.genes_pool[i]
    return child

def repopulate(population, parents, children_count):
    """Repopulate the shorten population by crossing the parents genes and adding new formed children""" 
    population_size = len(population)
    while children_count<population_size:
        p1_pos = get_parent_index(parents,None)
        p2_pos = get_parent_index(parents,p1_pos)
        population[children_count]=cross(parents[p1_pos],parents[p2_pos])
        population[children_count+1]=cross(parents[p2_pos],parents[p1_pos])
        children_count+=2


def mutate(population, chromosome_count, gene_count):
    """Mutate some of the genes of the specific amount of the chromosomes of given population""" 
    pop_size = len(population)
    for i in range(chromosome_count):
        chromosome_pos = random.randint(0, pop_size - 1)
        chromo = population[chromosome_pos]
        for j in range(gene_count): #Just selecting gene_count of genes and generating new random chord for them
            rand_note  = random.randint(0,max_note)
            rand_chord = Chord(rand_note, random.choice(CHORD_ALL))
            gene_pos = random.randint(0, chromo.size - 1)
            chromo.genes_pool[gene_pos] = rand_chord

iteration_count=0
max_iteration_count=5000
""" The main art of geneting algorithm where we create population, select the best of them, repopulate, mutate it in loop until we get the best result. """
while True: 
    iteration_count += 1
    calculate_rating(population)
    population = sorted(population)
    if(population[0].rating==0 or iteration_count>max_iteration_count):
        break
    select(population,survivors)
    size=len(population)
    repopulate(population,survivors, population_size//2)
    mutate(population,size//2,1) #We mutate always only one gene, as it will allow quicker evolution for almost perfect chords

"""Now we append track using the best generated accompanement and write it into the new output file"""
for chord in population[0].genes_pool:
      tracks.append(Message('note_on',channel=0,note=chord.note_list[0]+average_displacement,velocity=velocity,time=0))
      tracks.append(Message('note_on',channel=0,note=chord.note_list[1]+average_displacement,velocity=velocity,time=0))
      tracks.append(Message('note_on',channel=0,note=chord.note_list[2]+average_displacement,velocity=velocity,time=0))
      tracks.append(Message('note_off',channel=0,note=chord.note_list[0]+average_displacement,velocity=velocity,time=CHORD_DURATION))
      tracks.append(Message('note_off',channel=0,note=chord.note_list[1]+average_displacement,velocity=velocity,time=0))
      tracks.append(Message('note_off',channel=0,note=chord.note_list[2]+average_displacement,velocity=velocity,time=0))

midi.tracks.append(tracks)
midi.save(OUTPUT_FILE_NAME)