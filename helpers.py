from pathlib import Path
import pandas as pd
dif_percent = 30

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
    def get_from_te(fasta, te_name):
        te_len = Helpers.get_len_fasta(fasta)[te_name]
        return int(te_len - te_len * dif_percent / 100)

    @staticmethod
    def get_to_te(fasta, te_name):
        te_len = Helpers.get_len_fasta(fasta)[te_name]
        return int(te_len + te_len * dif_percent / 100)

    @staticmethod
    def read_bl(filename):
        header_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'expect value', 'bitscore']
        return pd.read_csv(filename, sep="\t", header=None, names=header_names)

    @staticmethod
    def save_df(df, filename):
        df.to_csv(filename, sep='\t', index=False, encoding='utf-8', header=False)

    @staticmethod
    def get_filebase(file):
        return Path(file).stem
