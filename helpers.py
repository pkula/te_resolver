from pathlib import Path
import pandas as pd
from const import BLAST_HEADER_NAMES

class Helpers:
    @staticmethod
    def get_fasta_files(directory):
        return [a for a in Path(directory).glob("*.fasta")] + [
            a for a in Path(directory).glob("*.fa")
        ]
    @staticmethod
    def get_fasta(filename):
        n = 0
        sequence = ""
        name = ""
        with open(filename) as f:
            for line in f.readlines():
                if line.startswith('>'):
                    n = n + 1
                    if not n == 1:
                        yield name, sequence
                    name = line[1:].strip()
                    sequence = ""
                else:
                    if name:
                        sequence = sequence + line.strip()
        yield name, sequence

    @staticmethod
    def get_len_fasta(filename):
        counter = {}
        with open(filename) as f:
            name = ""
            for line in f.readlines():
                if line.startswith('>'):
                    name = line[1:].strip()
                    counter[name] = 0
                else:
                    if name:
                        counter[name] = counter[name] + len(line.strip()) - line.count("~")
        return counter

    @staticmethod
    def get_from_te(config):
        te_len = Helpers.get_len_fasta(config.te_filepath)[config.te_name]
        return int(te_len - te_len * config.dif_percent / 100)

    @staticmethod
    def get_to_te(config):
        te_len = Helpers.get_len_fasta(config.te_filepath)[config.te_name]
        return int(te_len + te_len * config.dif_percent / 100)

    @staticmethod
    def read_bl(filename):
        return pd.read_csv(filename, sep="\t", header=None, names=BLAST_HEADER_NAMES)

    @staticmethod
    def save_df(df, filename):
        df.to_csv(filename, sep='\t', index=False, encoding='utf-8', header=False)

    @staticmethod
    def get_filebase(file):
        return Path(file).stem
