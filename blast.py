import logging
import os
from pathlib import Path
import pandas as pd
from helpers import Helpers
from config import Config


class Blast:
    def __init__(self):
        self.config = Config()

    # todo
    def make_ont_db(self):
        logging.info("start make db")
        pass_path = self.config.pass_reads_merged_path

        for f in Helpers.get_fasta_files(pass_path):
            name = pass_path / f.name.split(".")[0]
            os.system(f"makeblastdb -in {f} -dbtype nucl -out {name}")
        logging.info("end make db")

    #new
    @staticmethod
    def make_db(fasta_file, out):
        os.system(f"makeblastdb -in {fasta_file} -dbtype nucl -out {out}")

    @staticmethod
    def make_ont_db(config):
        for base, file in zip(config.ont_bases, config.ont_files):
            Blast.make_db(file, config.ont_path / base)

    #new
    @staticmethod
    def _run_te_ont(query, db, out):
        os.system(f"blastn -num_threads 20 -outfmt 6 -query {query} -db {db} -out {out} -dust no -perc_identity 0.9")
        return Helpers.read_bl(out)

    @staticmethod
    def run_te_ont(config):
        for base in config.ont_bases:
            Blast._run_te_ont(
                config.te_filepath,
                config.ont_path / base,
                config.get_te_ont_bl_path(base),
                )

    @staticmethod
    def run(query, db, out):
        os.system(f"blastn -num_threads 20 -outfmt 6 -query {query} -db {db} -out {out} -dust no -perc_identity 0.9")
        return Helpers.read_bl(out)
