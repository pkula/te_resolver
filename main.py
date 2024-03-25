import argparse
import logging
import os
from pathlib import Path

import pandas as pd

MAIN_PATH = "."

TE_BLAST_IDENT_PRC = 90
GENOME_BLAST_IDENT_PRC = 80

TE_IN = ""
ONT_IN = ""
GENOME_IN = ""

class Blast:
    def __init__(self):
        self.config = Config()

    @staticmethod
    def _get_fasta_files(directory):
        return [a for a in Path(directory).glob("*.fasta")] + [
            a for a in Path(directory).glob("*.fa")
        ]

    def make_pass_db(self):
        logging.info("start make db")
        pass_path = self.config.pass_reads_merged_path

        for f in self._get_fasta_files(pass_path):
            name = pass_path / f.name.split(".")[0]
            os.system(f"makeblastdb -in {f} -dbtype nucl -out {name}")
        logging.info("end make db")

    # todo not use
    def _run_blast(self, path, query, db, out):
        cmd = (
            f"blastn -num_threads 20 -outfmt 6 -query {query} "
            f"-db {db} -out {out} -dust no -perc_identity {TE_BLAST_IDENT_PRC/100}"
        )
        os.system(f"cd {path} && {cmd}")

    def run_pass_blast(self):
        logging.info("run blast")

        for f in self._get_fasta_files(self.config.pass_reads_merged_path):
            name = f.name.split(".")[0]
            name_path = self.config.pass_reads_merged_path / name
            te_ont_path = self.config.blast_path / f"TE_{name}.bl"
            # todo delete Alex1_Alex3.fasta
            cmd = (
                f"blastn -num_threads 20 -outfmt 6 -query /home/kira/projects/src/results/other_data/Alex1_Alex3.fasta "
                f"-db {name_path} -out {te_ont_path} -dust no -perc_identity {TE_BLAST_IDENT_PRC/100}"
            )
            os.system(cmd)
        logging.info("end run blast")

    # todo
    def run_genome_blast(self):
        cmd = f"makeblastdb -in {self.config.main_path}/genome/DCARv3.4.fa  -dbtype nucl -out {self.config.main_path}/genome/DCARv3.4"
        os.system(cmd)
        cmd = f"blastn -num_threads 20 -outfmt 6 -query {self.config.pass_reads_merged_path}/TE_K10B.masked.fasta -db {self.config.genome_path}/DCARv3.4 -out {self.config.blast_path}/TE_K10B_genome.bl -perc_identity {TE_BLAST_IDENT_PRC/100} -max_target_seqs 4"
        os.system(cmd)


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

        # set needed paths
        self.main_path = self.path = Path(MAIN_PATH).absolute() / "src"

        self.results_path = self.main_path / "results"
        self.other_data_path = self.results_path / "other_data"
        self.pass_reads_merged_path = self.results_path / "pass_reads_merged"

        self.analysis_path = self.main_path / "analysis"
        self.blast_path = self.analysis_path / "blast"

        self.genome_path = self.main_path / "genome"

        # Long reads from ONT sequencing files/directories -o
        self.raw_ont_path = Path("/home/kira/data2/ont")
        # Transcposon sequences (TE) files/directories -t
        self.raw_te_path = Path("/home/kira/data2/te")
        # Reference genome files/directories -g
        self.raw_genome_path = Path("/home/kira/data2/genome")

        self.create_directories()
        logging.info("end config")

    def create_directories(self):
        self.main_path.mkdir(parents=True, exist_ok=True)
        self.results_path.mkdir(parents=True, exist_ok=True)
        self.other_data_path.mkdir(parents=True, exist_ok=True)
        self.pass_reads_merged_path.mkdir(parents=True, exist_ok=True)
        self.analysis_path.mkdir(parents=True, exist_ok=True)
        self.blast_path.mkdir(parents=True, exist_ok=True)
        self.genome_path.mkdir(parents=True, exist_ok=True)


class Runner:
    def __init__(self):
        self.config = Config()

    @staticmethod
    def from_fq_to_fa(fq):
        fq = Path(fq).absolute()
        basename = Path(fq).name.split(".")[0]
        os.system(f"seqtk seq -a {fq} > {fq.parent / basename}.fasta")

    @staticmethod
    def _copy_file(from_file, to_file):
        os.system(f"cp {from_file} {to_file}")

    def copy_files(self):
        """
        # for all sequencing results
        #cat  FAS82558_pass*.fastq.gz > LK5_1_3.fastq.gz
        """
        # copy fastq.gz files into pass_reads_merged
        logging.info("start copying fq ont")
        for f in self.config.raw_ont_path.glob("*.fastq.gz"):
            self._copy_file(f, self.config.pass_reads_merged_path)
        logging.info("start copying fq ont")

        # upload Alex1_Alex3.fasta into other_data
        logging.info("start copying te")
        for f in self.config.raw_te_path.glob("*.fasta"):
            self._copy_file(f, self.config.other_data_path)
        logging.info("end copying te")

        logging.info("start copying genome")
        for f in self.config.raw_genome_path.glob("*"):
            self._copy_file(f, self.config.genome_path)
        logging.info("end copying te")

    def create_bed_files(self):
        """
        for file in *.bl; do
        awk 'OFS="\t"  {
        if ($9 <$10)
            print $2,$9,$10,$1;
        else
            print $2,$10,$9,$1;
        }' $file | sort -k1,1 -k2,2n |bedtools merge -i - -d 1500 -c 4 -o collapse  >$file.bed;
        done
        """

        def sort_row(row):
            if row["sstart"] > row["send"]:
                row["sstart"], row["send"] = row["send"], row["sstart"]
            return row

        column_names = [
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
        ]

        blast_path = self.config.blast_path
        for f in blast_path.glob("*.bl"):
            basename = Path(f.name).stem
            merged_path = str(blast_path / f"{basename}.merged")
            bed_path = str(blast_path / f"{f}.bed")

            csv = pd.read_csv(f"{f}", sep="\t", header=None, names=column_names)
            csv.apply(sort_row, axis=1).to_csv(
                merged_path,
                sep="\t",
                header=None,
                index=False,
                columns=["sseqid", "sstart", "send", "qseqid"],
            )
            os.system(
                f"cat {merged_path} | sort -k1,1 -k2,2n | bedtools merge -i - -d 1500 -c 4 -o collapse  > {bed_path}"
            )

    # todo
    def create_seq_list(self):
        cmd = f"for file in {self.config.blast_path}/*.bed; do awk '{{print $1}}' $file|sort|uniq > $file.seq_list.txt; done "
        os.system(cmd)
        cmd = f"for file in {self.config.blast_path}/*.txt; do wc -l $file >> {self.config.blast_path}/count.seq_list; done "
        os.system(cmd)

    # todo
    def subseq(self):
        for f in self.config.pass_reads_merged_path.glob("*.fasta"):
            name = f.name.split(".")[0]
            cmd = f"seqtk subseq {f} {self.config.blast_path}/TE_{name}.bl.bed.seq_list.txt >{self.config.pass_reads_merged_path}/TE_{name}.fasta"
            os.system(cmd)

    # todo
    def mask(self):
        cmd = f"bedtools maskfasta -fi {self.config.pass_reads_merged_path}/TE_K10B.fasta -bed {self.config.analysis_path}/blast/TE_K10B.bl.bed -fo {self.config.pass_reads_merged_path}/TE_K10B.masked.fasta"
        os.system(cmd)

    def convert_pass(self):
        logging.info("start converting pass from fastaq to fasta")
        """
        # change fq.gz into fasta (for blast)
        #seqtk you can install by: https://anaconda.org/bioconda/seqtk
        seqtk seq -a K10B.fastq.gz >K10B.fasta
        grep '>' K10B.fasta|wc -l #to display number of seq     6939343
        """
        pass_path = self.config.pass_reads_merged_path
        for f in pass_path.glob("*.fastq.gz"):
            self.from_fq_to_fa(f)
        logging.info("end converting pass to fasta")


if __name__ == "__main__":
    runner = Runner()
    blast = Blast()
    runner.copy_files()
    runner.convert_pass()
    blast.make_pass_db()
    blast.run_pass_blast()
    runner.create_bed_files()
    runner.create_seq_list()
    runner.subseq()
    runner.mask()
    blast.run_genome_blast()

    #os.system(f"rm -rf {runner.config.main_path}")
