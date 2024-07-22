import random

from src.helpers import Helpers
from bioinfokit.analys import Fasta


def _rand_split(sequence, splits):
    for splitLen in splits:
        if splitLen > len(sequence):
            break
        yield sequence[:splitLen]
        sequence = sequence[splitLen::]


def _rand_gen(int_list):
    while True:
        yield int(random.choice(int_list))


def cut_genome(genome_path, ont_path, coverage, filename, from_n):
    possibility_len = list(Helpers.get_len_fasta(ont_path).values())
    n = 0
    for i in range(from_n, from_n + coverage):
        print(i)
        for name, genome_sequence in Fasta.fasta_reader(file=genome_path):
            print(name)
            with open(filename, "a") as file:
                for seq in _rand_split(genome_sequence, _rand_gen(possibility_len)):
                    n = n + 1
                    file.write(f">{name}_{i}_{n}\n{seq}\n")
