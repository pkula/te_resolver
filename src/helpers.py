from pathlib import Path

import pandas as pd
from bioinfokit.analys import Fasta

from src.const import BLAST_HEADER_NAMES


class Helpers:
    @staticmethod
    def get_fasta_files(directory):
        return [a for a in Path(directory).glob("*.fasta")] + [
            a for a in Path(directory).glob("*.fa")
        ]

    @staticmethod
    def get_len_fasta(filename):
        return {
            header: len(sequence)
            for header, sequence in Fasta.fasta_reader(file=filename)
        }

    @staticmethod
    def get_from_te(config):
        return int(config.te_len - config.te_len * config.dif_percent / 100)

    @staticmethod
    def get_to_te(config):
        return int(config.te_len + config.te_len * config.dif_percent / 100)

    @staticmethod
    def read_bl(filename):
        return pd.read_csv(filename, sep="\t", header=None, names=BLAST_HEADER_NAMES)

    @staticmethod
    def save_df(df, filename):
        df.to_csv(filename, sep="\t", index=False, encoding="utf-8", header=False)

    @staticmethod
    def get_filebase(file):
        return Path(file).stem
