import argparse
import os
import logging
import pandas as pd
from pathlib import Path


class Blast:
    def __init__(self):
        self.config = Config()

    # blast for sequence Alex1_Alex3

    def make_db(self):
        logging.info("start make db")
        # 1) make a database from your nucleotide sequences.
        pass_path = self.config.pass_reads_merged_path
        os.system(
            f"makeblastdb -in {pass_path}/LK5_1_3.fasta -dbtype nucl -out {pass_path}/LLK5_1_3"
        )
        os.system(
            f"makeblastdb -in {pass_path}/K10B.fasta  -dbtype nucl -out {pass_path}/K10B"
        )
        os.system(
            f"makeblastdb -in {pass_path}/K10f.fasta  -dbtype nucl -out {pass_path}/K10f"
        )
        logging.info("end make db")

    def run_blast(self):
        logging.info("run blast")

        main_path = self.config.main_path
        analysis_path = self.config.analysis_path

        # 2) run blast https://www.ncbi.nlm.nih.gov/books/NBK279684/) from /home/MSA1/username
        # todo k value
        cmd = (
            "blastn -num_threads 20 -outfmt 6 -query results/other_data/Alex1_Alex3.fasta "
            "-db results/pass_reads_merged/LK5_1_3 -out analysis/blast/Alex_LK5_1_3.bl -dust no -perc_identity 0.9"
        )
        os.system(f"cd {main_path} && {cmd}")
        cmd = (
            "blastn -num_threads 20 -outfmt 6 -query results/other_data/Alex1_Alex3.fasta "
            "-db results/pass_reads_merged/K10B -out analysis/blast/Alex_K10B.bl -dust no -perc_identity 0.9"
        )
        os.system(f"cd {main_path} && {cmd}")
        cmd = (
            "blastn -num_threads 20 -outfmt 6 -query results/other_data/Alex1_Alex3.fasta "
            "-db results/pass_reads_merged/K10f -out analysis/blast/Alex_K10f.bl -dust no -perc_identity 0.9"
        )
        os.system(f"cd {main_path} && {cmd}")
        logging.info("end run blast")


class Config:
    def __new__(cls):
        if not hasattr(cls, "instance"):
            cls.instance = super(Config, cls).__new__(cls)
        return cls.instance

    def __init__(self):
        logging.info("start config")
        # prepare parser
        parser = argparse.ArgumentParser(description="TE pipeline resolver")

        parser.add_argument("-w", action="store_true", help="flag default false")
        parser.add_argument("-k", action="store", type=int, default=90, help="")
        parser.add_argument("-v", "--verbose", action="store_true")

        self.path = Path(".").absolute()  # /home/username/projects/te_resolver

        # set needed paths

        # /home/username/projects
        self.main_path = self.path.parent  # m main path -> analysis and results
        self.results_path = (
            self.main_path / "results"
        )  # /home/username/projects/results
        self.other_data_path = (
            self.results_path / "other_data"
        )  # /home/username/projects/results/other_data
        self.pass_reads_merged_path = (
            self.results_path / "pass_reads_merged"
        )  # /home/username/projects/results/pass_reads_merged

        self.analysis_path = (
            self.main_path / "analysis"
        )  # /home/username/projects/analysis
        self.blast_path = (
            self.analysis_path / "blast"
        )  # /home/username/projects/analysis




        # create needed paths



        # mkdir m paths
        self.main_path.mkdir(parents=True, exist_ok=True)
        self.results_path.mkdir(parents=True, exist_ok=True)
        self.other_data_path.mkdir(parents=True, exist_ok=True)
        self.pass_reads_merged_path.mkdir(parents=True, exist_ok=True)
        self.analysis_path.mkdir(parents=True, exist_ok=True)
        self.blast_path.mkdir(parents=True, exist_ok=True)






        self.args = parser.parse_args()

        # Long reads from ONT sequencing files/directories -o
        self.ont_path = Path("/home/kira/data2/ont")
        # Transcposon sequences (TE) files/directories -t
        self.te_path = Path("/home/kira/data2/te")
        # Reference genome files/directories -g
        self.genome_path = Path("/home/kira/data2/genome")

        logging.info("end config")


class Runner:
    def __init__(self):
        self.config = Config()

    @staticmethod
    def from_fq_to_fa(fq, fa):
        os.system(f"seqtk seq -a {fq} > {fa}")

    @staticmethod
    def _copy_file(from_file, to_file):
        os.system(f"cp {from_file} {to_file}")

    def preparation(self):
        logging.info("start preparation")
        """
        # for all sequencing results
        #cat  FAS82558_pass*.fastq.gz > LK5_1_3.fastq.gz
        """
        other_data_path = self.config.other_data_path
        pass_reads_merged_path = self.config.pass_reads_merged_path
        # copy fastq.gz files into pass_reads_merged
        for f in self.config.ont_path.glob("*.fastq.gz"):
            self._copy_file(f, pass_reads_merged_path)
        # upload Alex1_Alex3.fasta into other_data
        for f in self.config.genome_path.glob("*.fasta"):
            self._copy_file(f, other_data_path)
        logging.info("end preparation")

    def create_bed_files(self):
        """
        awk 'OFS="\t"  {
        if ($9 <$10)
             print $2,$9,$10,$1;
        else
            print $2,$10,$9,$1;
        }' $file | sort -k1,1 -k2,2n |bedtools merge -i - -d 1500 -c 4 -o collapse  >$file.bed;
        done
        """

        def xxx(row):
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
            
            okno = pd.read_csv(
                f"{f}", sep="\t", header=None, names=column_names
            )
            bed_path = str(blast_path / f"{f.name}.bed")
            bed_path2 = str(blast_path / f"{f.name}.bed2")
            okno.apply(xxx, axis=1).to_csv(
                bed_path,
                sep="\t",
                header=None,
                index=False,
                columns=["sseqid", "sstart", "send", "qseqid"],
            )
            os.system(
                f"sort -k1,1 -k2,2n {bed_path} >{bed_path2}"
            )
            merged_path = str(blast_path / f"{f.name}.merged")
            cmd = f"bedtools merge -i {bed_path2} -d 1500 -c 4 -o collapse >{merged_path}"
            os.system(
                f"bedtools merge -i {bed_path2} -d 1500 -c 4 -o collapse >{merged_path}"
            )


    def preparation2(self):
        logging.info("start preparation2")
        """
        #fq.gz into fasta (for blast)
        #seqtk you can install by: https://anaconda.org/bioconda/seqtk
        seqtk seq -a LK5_1_3.fastq.gz > LK5_1_3.fasta
        grep '>' LK5_1_3.fasta|wc -l #to display number of seq  1393171
        seqtk seq -a K10B.fastq.gz >K10B.fasta
        grep '>' K10B.fasta|wc -l #to display number of seq     6939343
        seqtk seq -a K10f.fastq.gz >K10f.fasta
        grep '>' K10f.fasta|wc -l #to display number of seq     4775182
        """
        pass_path = self.config.pass_reads_merged_path
        fq1 = pass_path / "LK5_1_3.fastq.gz"
        fa1 = pass_path / "LK5_1_3.fasta"
        self.from_fq_to_fa(fq1, fa1)
        fq2 = pass_path / "K10B.fastq.gz"
        fa2 = pass_path / "K10B.fasta"
        self.from_fq_to_fa(fq2, fa2)
        fq3 = pass_path / "K10f.fastq.gz"
        fa3 = pass_path / "K10f.fasta"
        self.from_fq_to_fa(fq3, fa3)

        logging.info("end preparation2")


if __name__ == "__main__":
    runner = Runner()
    runner.preparation()
    runner.preparation2()
    blast = Blast()
    blast.make_db()
    blast.run_blast()

    #runner.create_bed_files()
    import ipdb; ipdb.set_trace()

    os.system(f"rm -rf {results_path}")
    os.system(f"rm -rf {analysis_path}")
