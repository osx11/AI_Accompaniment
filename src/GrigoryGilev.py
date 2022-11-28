import random
from typing import List, Union
import mido
from mido import Message, MetaMessage
import music21

INPUT_FILE_NAME = 'input1.mid'
OUTPUT_FILE_NAME = 'output1.mid'

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


def avg_velocity() -> int:
    """
    :returns: avg sound level of original song
    :rtype: int
    """

    notes_list = [note.velocity for note in midi.tracks[1] if isinstance(note, Message) and note.type == 'note_on']
    return int(sum(notes_list)/len(notes_list))


def avg_octave() -> int:
    """
    :returns: avg notes range of original sound
    :rtype: int
    """

    notes_list = [int(note.note/12) for note in midi.tracks[1] if isinstance(note, Message) and note.type == 'note_on']
    return int(sum(notes_list)/len(notes_list))


class Chord:
    def __init__(self, note_start: int, chord_type: List[int]):
        """
        :param note_start: start note
        :param chord_type: array of notes (representing its type (MAJOR, MINOR, DIM))
        """

        self.root_note = note_start
        self.type = chord_type
        self.note_list = [(note_start + note) % 12 for note in chord_type]

    def has_note(self, note: int) -> bool:
        """
        :param note: note to check
        :returns: true if note is present in the current chord, otherwise false
        :rtype: bool
        """
        if note is None:
            return False

        return any(chord_note == note % 12 for chord_note in self.note_list)

    def __eq__(self, other):
        return self.root_note == other.root_note and self.type == other.type


class Accompaniment:
    def __init__(self, scale: int, ticks: int):
        """
        :param scale: scale of the song
        :param ticks: number of ticks (tacts) to calculate the nubmer of chords
        """

        self.scale = scale
        self.ticks = ticks
        self.accordantChords = []
        self.accordantChords.append(Chord(scale % 12, CHORD_MAJOR))
        self.accordantChords.append(Chord((scale + 5) % 12, CHORD_MAJOR))
        self.accordantChords.append(Chord((scale + 7) % 12, CHORD_MAJOR))
        self.accordantChords.append(Chord((scale + 9) % 12, CHORD_MINOR))
        self.accordantChords.append(Chord((scale + 2) % 12, CHORD_MINOR))
        self.accordantChords.append(Chord((scale + 4) % 12, CHORD_MINOR))
        self.accordantChords.append(Chord((scale + 11) % 12, CHORD_DIM))

    def get_accordant_chord(self, note: int) -> Chord:
        """
        :param note: note to find the fitting chord
        :returns: chord with the note
        :rtype: Chord
        """

        for _chord in self.accordantChords:
            if _chord.has_note(note):
                return _chord

        return random.choice(self.accordantChords)

    def has_in_accordant_chords(self, _chord: Union[Chord, None]) -> bool:
        """
        :param _chord: chord to check
        :returns: true of the accordant chord contains the note, otherwise false
        :rtype: bool
        """

        return any(existing_chord == _chord for existing_chord in self.accordantChords)

    """ Return if the accordant chords contain the given note inside any of them"""
    def has_note_in_accordant_chords(self, note: Union[int, None]):
        """
        :param note: note to check
        :returns: true if accordant chords has the note, otherwise false
        :rtype: bool
        """

        if note is None:
            return False

        return any(existing_chord.has_note(note % 12) for existing_chord in self.accordantChords)


def get_notes_amount(_track) -> int:
    """
    :returns: number of half-tracks of the track
    :rtype: int
    """

    beats = sum(_msg.time for _msg in _track if type(_msg) is Message)
    length = (beats + CHORD_DURATION - 1) // CHORD_DURATION

    return length


def compute_beats(_track) -> List[int]:
    """
    :returns: notes of the track
    :rtype:
    """

    length = get_notes_amount(_track)
    notes_list = [0]*length
    beats = 0
    last_note = 0

    for _msg in _track:
        if type(_msg) is Message:
            if beats % CHORD_DURATION == 0 and _msg.type == "note_on":
                if notes_list[beats//CHORD_DURATION] is None:
                    notes_list[beats//CHORD_DURATION] = _msg.note

            if _msg.type == "note_off":
                last_note = _msg.note

            beats += _msg.time

    if beats % CHORD_DURATION == 0:
        notes_list[-1] = last_note

    return notes_list


tempo = [msg.tempo for track in midi.tracks for msg in track if isinstance(msg, MetaMessage) and msg.type == "set_tempo"][-1]
tracks = [MetaMessage('set_tempo', tempo=tempo, time=0)]

song_notes = compute_beats(midi.tracks[1])
velocity = int(avg_velocity() * 0.9)
avg_displacement = 12 * (avg_octave() - 1)

final_genome_chord = Accompaniment(MUSIC_SCALE, len(song_notes))
MAX_NOTE = 120


class Chromosome:
    def __init__(self, ticks):
        """
        :param ticks: number of genera in accompaniment
        """

        self.ticks = ticks
        self.genes_pool = [None]*ticks
        self.rating = 0
        self.generate_random_genes()

    def generate_random_genes(self):
        """
        Random genes generator. Uses note in range [0, 120]
        """

        for i in range(self.ticks):
            rand_note = random.randint(0, MAX_NOTE)
            rand_chord = Chord(rand_note, random.choice(CHORD_ALL))
            self.genes_pool[i] = rand_chord

    def __eq__(self, other):
        return self.rating == other.rating

    def __lt__(self, other):
        return self.rating < other.rating


def create_population(population_size: int, chromosome_size: int) -> List[Union[None, Chromosome]]:
    """
    :param population_size: size of population
    :param chromosome_size: number of chromosomes
    :returns: chromosomes population. genes amount is fixed
    :rtype: List[Union[None, Chromosome]]
    """

    _population = [None]*population_size

    for i in range(population_size):
        _population[i] = Chromosome(chromosome_size)

    return _population


def calculate_rating(_population: List[Union[None, Chromosome]]):
    """
    Fitness function to calculate ranting of chromosomes
    """
    for chromosome in _population:
        chromosome.rating = chromosome.ticks
        for i in range(chromosome.ticks):
            if final_genome_chord.has_in_accordant_chords(chromosome.genes_pool[i]):
                chromosome.rating -= 0.5
                if song_notes[i] is None:
                    chromosome.rating -= 0.5
                    continue
                if not final_genome_chord.has_note_in_accordant_chords(song_notes[i]) or chromosome.genes_pool[i].has_note(song_notes[i]):
                    chromosome.rating -= 0.5


POPULATION_SIZE = 128
best_fit_chromosomes = [None] * (POPULATION_SIZE // 4)

population = create_population(POPULATION_SIZE, final_genome_chord.ticks)


def select(_population: List[Union[None, Chromosome]], survivors):
    """
    Selects top of best fit chromosomes
    """

    for i in range(len(survivors)):
        survivors[i] = _population[i]


def get_parent_index(parents: List[Chromosome], other_parent_index: Union[int, None]) -> int:
    """
    :returns: index of random chromosome from parent such that the other parent is not that index
    """

    while True:
        index = random.randint(0, len(parents)-1)
        if other_parent_index is None or other_parent_index != index:
            return index


def cross(chromosome1: Chromosome, chromosome2: Chromosome) -> Chromosome:
    """
    Crossover algorithm
    :returns: child of two chromosomes
    :rtype: Chromosome
    """

    ticks = chromosome1.ticks
    point = random.randint(0, ticks-1)
    child = Chromosome(ticks)

    for i in range(point):
        child.genes_pool[i] = chromosome1.genes_pool[i]

    for i in range(point, ticks):
        child.genes_pool[i] = chromosome2.genes_pool[i]

    return child


def populate(_population, parents: List[Union[Chromosome, None]], children_count: int):
    population_size = len(_population)
    while children_count < population_size:
        p1_pos = get_parent_index(parents, None)
        p2_pos = get_parent_index(parents, p1_pos)
        _population[children_count] = cross(parents[p1_pos], parents[p2_pos])
        _population[children_count+1] = cross(parents[p2_pos], parents[p1_pos])
        children_count += 2


def mutate(_population: List[Union[None, Chromosome]], chromosome_count: int, gene_count: int):
    pop_size = len(_population)
    for i in range(chromosome_count):
        chromosome_pos = random.randint(0, pop_size - 1)
        chromosome = _population[chromosome_pos]
        for j in range(gene_count):
            rand_note = random.randint(0, MAX_NOTE)
            rand_chord = Chord(rand_note, random.choice(CHORD_ALL))
            gene_pos = random.randint(0, chromosome.ticks - 1)
            chromosome.genes_pool[gene_pos] = rand_chord


iteration_count = 0
max_iteration_count = 5000

while True:
    iteration_count += 1
    calculate_rating(population)
    population = sorted(population)

    if population[0].rating == 0 or iteration_count>max_iteration_count:
        break

    select(population, best_fit_chromosomes)
    size = len(population)
    populate(population, best_fit_chromosomes, POPULATION_SIZE // 2)
    mutate(population, size // 2, 1)

"""Now we append track using the best generated accompanement and write it into the new output file"""
for chord in population[0].genes_pool:
    tracks.append(Message('note_on', channel=0, note=chord.note_list[0] + avg_displacement, velocity=velocity, time=0))
    tracks.append(Message('note_on', channel=0, note=chord.note_list[1] + avg_displacement, velocity=velocity, time=0))
    tracks.append(Message('note_on', channel=0, note=chord.note_list[2] + avg_displacement, velocity=velocity, time=0))
    tracks.append(Message('note_off', channel=0, note=chord.note_list[0] + avg_displacement, velocity=velocity, time=CHORD_DURATION))
    tracks.append(Message('note_off', channel=0, note=chord.note_list[1] + avg_displacement, velocity=velocity, time=0))
    tracks.append(Message('note_off', channel=0, note=chord.note_list[2] + avg_displacement, velocity=velocity, time=0))

midi.tracks.append(tracks)
midi.save(OUTPUT_FILE_NAME)
