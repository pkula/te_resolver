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

        # main filders
        self.main_path = Path(MAIN_PATH).absolute()
        self.results_path = self.main_path / "results"

        # input data
        self.genome_path = self.main_path / "genome"
        self.ont_path = self.main_path / "ont"
        self.te_path = self.main_path / "te"

        # args
        self.te_name = TE_NAME
        self.te_filepath = self.te_path / TE_FILE
        self.genome_filepath = self.genome_path / GENOME_FILE

        # first part
        self.first_path = self.results_path / "first"
        self.first_blast_path = self.first_path / "blast"
        self.first_bed_path = self.first_path / "bed"
        self.first_subseq_path = self.first_path / "subseq"

        # second part
        self.second_path = self.results_path / "second"
        self.report_filebase = f"{self.genome_filepath.stem}_{'_'.join(self.ont_bases)}"
        self.filtered_records_filepath = self.second_path / f"filtered_{self.report_filebase}"
        self.raw_report_filepath = self.second_path / f"raw_report_{self.report_filebase}"
        self.final_report_filepath = self.second_path / f"final_report_{self.report_filebase}"

        # masked genome
        self.masked_genome_path = self.results_path / "masked_genome"

        # set needed paths
        self.te_len = Helpers.get_len_fasta(self.te_filepath)[self.te_name]
        self.from_te = Helpers.get_from_te(self.te_filepath, self.te_name)
        self.to_te = Helpers.get_to_te(self.te_filepath, self.te_name)

        print(f"ONT: {self.ont_bases}")
        print(f"Genome: {self.genome_filepath}")
        print(f"TE: {self.te_name}")

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
