import argparse
import logging
from pathlib import Path

from bioinfokit.analys import Fasta

from src.const import (BLAST_THREADS, DIF_PERCENT, FLANKS_LEN, GROUP_LEN,
                       MAIN_PATH, MAX_NONREF_LEN)
from src.helpers import Helpers


class Config:
    def __new__(cls):
        if not hasattr(cls, "instance"):
            cls.instance = super(Config, cls).__new__(cls)
        return cls.instance

    def prepare_parser(self):
        parser = argparse.ArgumentParser(description="TE pipeline resolver")
        parser.add_argument(
            "-s",
            "--script",
            action="store",
            type=str,
            default="finder", # finder, cutter,
            help="Script",
        )
        parser.add_argument(
            "-p",
            "--main_path",
            action="store",
            type=str,
            default=MAIN_PATH,
            help="Main path where ont, genome and te directories exist",
        )
        parser.add_argument("-t", "--te_file", action="store", help="TE filename")
        parser.add_argument(
            "-g", "--genome_file", action="store", help="Genome filename"
        )
        parser.add_argument("-e", "--te_name", action="store", help="TE name")
        parser.add_argument(
            "-f",
            "--flanks",
            action="store",
            type=int,
            default=FLANKS_LEN,
            help="Length of flanks",
        )
        parser.add_argument(
            "-r",
            "--groups",
            action="store",
            type=int,
            default=GROUP_LEN,
            help="Length of groups",
        )
        parser.add_argument(
            "-m",
            "--max_nonref_len",
            action="store",
            type=int,
            default=MAX_NONREF_LEN,
            help="Max nonref length",
        )
        parser.add_argument(
            "-d",
            "--dif_percent",
            action="store",
            type=int,
            default=DIF_PERCENT,
            help="Set length of TE from to",
        )
        parser.add_argument(
            "-b",
            "--blast_threads",
            action="store",
            type=int,
            default=BLAST_THREADS,
            help="Number of blast threads",
        )
        # todo not used - in run
        parser.add_argument(
            "-c", "--convert_fq", action="store_true", help="Convert ont fq to fasta"
        )
        parser.add_argument(
            "-l", "--create_db", action="store_false", help="Make blast dbs"
        )
        parser.add_argument("-o", "--do_blast", action="store_false", help="Do blasts")
        return parser

    def __init__(self):
        logging.info("start config")
        parser = self.prepare_parser()
        args = parser.parse_args()

        self.script = args.script

        self.convert_fq = args.convert_fq
        self.create_db = args.create_db
        self.do_blast = args.do_blast


        self.main_path = (
            Path(args.main_path).absolute()
            if args.main_path
            else Path(MAIN_PATH).absolute()
        )

        self.flanks_len = args.flanks
        self.group_len = args.groups
        self.max_nonref_len = args.max_nonref_len
        self.dif_percent = args.dif_percent
        self.blast_threads = args.blast_threads

        # main filders
        self.results_path = self.main_path / "results"

        # input data
        self.genome_path = self.main_path / "genome"
        self.ont_path = self.main_path / "ont"
        self.te_path = self.main_path / "te"
        self.db_path = self.main_path / "db"

        # args
        self.te_filepath = (
            self.te_path / args.te_file
            if args.te_file
            else [a for a in Path(self.te_path).glob("*")][0]
        )
        self.te_name = (
            args.te_name
            if args.te_name
            else Fasta.fasta_reader(file=self.te_filepath).__next__()[0]
        )
        self.genome_filepath = (
            self.genome_path / args.genome_file
            if args.genome_file
            else [a for a in Path(self.genome_path).glob("*")][0]
        )

        # first part
        self.first_path = self.results_path / "first"
        self.first_blast_path = self.first_path / "blast"
        self.first_bed_path = self.first_path / "bed"
        self.first_subseq_path = self.first_path / "subseq"

        # second part
        self.second_path = self.results_path / "second"
        self.report_filebase = (
            f"{self.genome_filepath.stem}_{'_'.join(self.ont_bases)}_{self.te_name}"
        )
        self.filtered_records_filepath = (
            self.second_path / f"filtered_{self.report_filebase}"
        )
        self.raw_report_filepath = (
            self.second_path / f"raw_report_{self.report_filebase}"
        )
        self.final_report_filepath = (
            self.second_path / f"final_report_{self.report_filebase}"
        )

        # masked genome
        self.masked_genome_path = self.results_path / "masked_genome"

        # set needed paths
        self.te_len = Helpers.get_len_fasta(self.te_filepath)[self.te_name]
        self.from_te = Helpers.get_from_te(self)
        self.to_te = Helpers.get_to_te(self)

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
        self.db_path.mkdir(parents=True, exist_ok=True)

    @property
    def ont_files(self):
        return Helpers.get_fasta_files(self.ont_path)

    @property
    def ont_bases(self):
        return [f.stem for f in self.ont_files]

    def get_te_ont_bl_path(self, ont_filebase):
        return (
            self.first_blast_path
            / f"TE_{self.genome_filepath.stem}_{ont_filebase}_{self.te_name}.bl"
        )

    def get_ont_filepath_from_ont_base(self, ont_filebase):
        for file in self.ont_files:
            if file.stem == ont_filebase:
                return file
