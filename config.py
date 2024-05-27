import argparse
import logging
from helpers import Helpers
from pathlib import Path

MAIN_PATH = "."

TE_FILE = "Alex1_3.fas"
TE_NAME = 'Alex1'
GENOME_FILE = "genomic2.fna"


class Config:
    def __new__(cls):
        if not hasattr(cls, "instance"):
            cls.instance = super(Config, cls).__new__(cls)
        return cls.instance

    def prepare_parser(self):
        parser = argparse.ArgumentParser(description="TE pipeline resolver")

        parser.add_argument("-w", action="store_true", help="flag default false")
        parser.add_argument("-k", action="store", type=int, default=90, help="")
        parser.add_argument("-v", "--verbose", action="store_true")
        return parser

    def __init__(self):
        logging.info("start config")
        parser = self.prepare_parser()
        self.args = parser.parse_args()

        self.main_path = Path(MAIN_PATH).absolute()
        self.genome_path = self.main_path / "genome"
        self.second_path = self.main_path / "second"
        self.ont_path = self.main_path / "ont"
        self.te_path = self.main_path / "te"
        self.masked_genome_path = self.main_path / "masked_genome"

        self.te_name = TE_NAME
        self.te_filepath = self.te_path / TE_FILE
        self.genome_filepath = self.genome_path / GENOME_FILE

        self.report_filebase = f"{self.genome_filepath.stem}_{'_'.join(self.ont_bases)}"
        self.filtered_records_filepath = self.main_path / f"filtered_{self.report_filebase}"
        self.raw_report_filepath = self.main_path / f"raw_report_{self.report_filebase}"
        self.final_report_filepath = self.main_path / f"final_report_{self.report_filebase}"


        # set needed paths


        self.results_path = self.main_path / "results"

        self.first_path = self.results_path / "first"
        self.first_blast_path = self.first_path / "blast"
        self.first_bed_path = self.first_path / "bed"
        self.first_subseq_path = self.first_path / "subseq"

        self.te_len = Helpers.get_len_fasta(self.te_filepath)[self.te_name]
        self.from_te = Helpers.get_from_te(self.te_filepath, self.te_name)
        self.to_te = Helpers.get_to_te(self.te_filepath, self.te_name)



        self.create_directories()
        logging.info("end config")

    def create_directories(self):
        self.results_path.mkdir(parents=True, exist_ok=True)
        self.first_path.mkdir(parents=True, exist_ok=True)
        self.first_blast_path.mkdir(parents=True, exist_ok=True)
        self.first_bed_path.mkdir(parents=True, exist_ok=True)
        self.first_subseq_path.mkdir(parents=True, exist_ok=True)
        self.masked_genome_path.mkdir(parents=True, exist_ok=True)
        self.second_path.mkdir(parents=True, exist_ok=True)


    @property
    def ont_files(self):
        return Helpers.get_fasta_files(self.ont_path)


    @property
    def ont_bases(self):
        return [f.stem for f in self.ont_files]
