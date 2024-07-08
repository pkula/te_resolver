import logging
import os
from pathlib import Path

import pandas as pd

from src.config import Config
from src.helpers import Helpers


class Blast:
    def __init__(self):
        self.config = Config()

    # new
    @staticmethod
    def make_db(fasta_file, out):
        os.system(f"makeblastdb -in {fasta_file} -dbtype nucl -out {out}")

    @staticmethod
    def make_ont_db(config):
        for base in config.ont_bases:
            file = config.get_ont_filepath_from_ont_base(base)
            Blast.make_db(file, config.db_path / base)

    # new
    @staticmethod
    def _run_te_ont(query, db, out, num_threads):
        os.system(
            f"blastn -num_threads {num_threads} -outfmt 6 -query {query} -db {db} -out {out} -dust no -perc_identity 0.9"
        )
        return Helpers.read_bl(out)

    @staticmethod
    def run_te_ont(config):
        for base in config.ont_bases:
            Blast._run_te_ont(
                config.te_filepath,
                config.db_path / base,
                config.get_te_ont_bl_path(base),
                config.blast_threads,
            )

    @staticmethod
    def run(query, db, out, num_threads):
        os.system(
            f"blastn -num_threads {num_threads} -outfmt 6 -query {query} -db {db} -out {out} -dust no -perc_identity 0.9"
        )
        return Helpers.read_bl(out)
