import argparse
import logging
import os
from pathlib import Path
import pandas as pd
from bioinfokit.analys import Fasta
from Bio import SeqIO
import numpy as np
from blast import Blast
from helpers import Helpers
from second_service import SecondService
from mask_genome import mask_genome


TE_BLAST_IDENT_PRC = 90
GENOME_BLAST_IDENT_PRC = 80

ont = ["ccc"]
ont = [
    "K10B",
    "K10f",
    "LK5_1_3"
]
genome_path = "genome/DCARv3.4.fa"
genome_path = "genome/genomic2.fna"
te_path = "te/Alex1_3.fas"
dif_percent = 30

#te_name = 'Alex1'
#te_len = Helpers.get_len_fasta(te_path)[te_name]
#from_te = Helpers.get_from_te(te_path, te_name)
#to_te = Helpers.get_to_te(te_path, te_name)


class Runner:
    def __init__(self):
        self.config = Config()

    @staticmethod
    def from_fq_to_fa(fq):
        fq = Path(fq).absolute()
        basename = Path(fq).name.split(".")[0]
        os.system(f"seqtk seq -a {fq} > {fq.parent / basename}.fasta")

    def correct_first_part(self, fq=True, db=True, do_blast=True):
        os.system("mkdir -p blast")
        os.system("mkdir -p first")
        if fq:
            for file in ont:
                os.system(f"seqtk seq -a ont/{file}.fastq.gz > ont/{file}.fasta")
        if db:
            print('okno')
            for file in ont:
                print(f"makeblastdb -in ont/{file}.fasta -dbtype nucl -out ont/{file}")
                os.system(f"makeblastdb -in ont/{file}.fasta -dbtype nucl -out ont/{file}")
        if do_blast:
            print('okno2')
            for file in ont:
                print(
                    f"blastn -num_threads 20 -outfmt 6 -query {te_path} -db ont/{file} " \
                    f"-out blast/TE_{file}_{te_name}.bl -dust no -perc_identity 0.9"
                )
                os.system(
                    f"blastn -num_threads 20 -outfmt 6 -query {te_path} -db ont/{file} " \
                    f"-out blast/TE_{file}_{te_name}.bl -dust no -perc_identity 0.9"
                )

        for o in ont:
            te_ont_bl = Helpers.read_bl(f"blast/TE_{o}_{te_name}.bl")
            te_ont_bl = te_ont_bl[te_ont_bl.qseqid.str.contains(te_name)]
            plus_te_ont_bl = te_ont_bl[te_ont_bl.sstart < te_ont_bl.send]
            minus_te_ont_bl = te_ont_bl[te_ont_bl.sstart > te_ont_bl.send]

            plus_te_ont_bed = plus_te_ont_bl[['sseqid', 'sstart', 'send', 'qseqid']].sort_values(
                by=['sseqid', 'sstart', 'send'])
            minus_te_ont_bed = minus_te_ont_bl[['sseqid', 'send', 'sstart', 'qseqid']].sort_values(
                by=['sseqid', 'send', 'sstart'])

            Helpers.save_df(minus_te_ont_bed, f'first/m_raw_first_blast_{o}_{te_name}.bed')
            Helpers.save_df(plus_te_ont_bed, f'first/p_raw_first_blast_{o}_{te_name}.bed')

            os.system(
                f"bedtools merge -i first/p_raw_first_blast_{o}_{te_name}.bed -d 1500 -c 4 -o collapse  > first/p_first_blast_{o}_{te_name}.bed.merged")
            os.system(
                f"bedtools merge -i first/m_raw_first_blast_{o}_{te_name}.bed -d 1500 -c 4 -o collapse  > first/m_first_blast_{o}_{te_name}.bed.merged")

            print("first part")

            plus_merged = pd.read_csv(
                f"first/p_first_blast_{o}_{te_name}.bed.merged", sep="\t", header=None).rename(
                {0: 'sseqid', 1: 'sstart', 2: 'send', 3: 'qseqid'}, axis=1
            )
            minus_merged = pd.read_csv(
                f"first/m_first_blast_{o}_{te_name}.bed.merged", sep="\t", header=None).rename(
                {0: 'sseqid', 1: 'sstart', 2: 'send', 3: 'qseqid'}, axis=1
            )

            plus_filtered_bed = plus_merged[(plus_merged.send - plus_merged.sstart) > from_te]
            minus_filtered_bed = minus_merged[(minus_merged.send - minus_merged.sstart) > from_te]

            # todo in this place

            Helpers.save_df(pd.DataFrame({'id': plus_filtered_bed['sseqid'].unique()}),
                               f'first/p_seq_list_{o}_{te_name}.txt')
            Helpers.save_df(pd.DataFrame({'id': minus_filtered_bed['sseqid'].unique()}),
                               f'first/m_seq_list_{o}_{te_name}.txt')
            # df.to_csv(filename, sep='\t', index=False, encoding='utf-8', header=False)
            print("second part")

            # long
            # it's filtered only our fasta
            os.system(
                f"seqtk subseq ont/{o}.fasta first/p_seq_list_{o}_{te_name}.txt > first/p_normal_{o}_{te_name}.fasta")
            os.system(
                f"seqtk subseq ont/{o}.fasta first/m_seq_list_{o}_{te_name}.txt > first/m_normal_{o}_{te_name}.fasta")
            print("1 part")

            create_subsequent_fasta(f"first/p_normal_{o}_{te_name}.fasta", plus_filtered_bed, f"p_{o}_{te_name}.fasta")
            create_subsequent_fasta(f"first/m_normal_{o}_{te_name}.fasta", minus_filtered_bed, f"m_{o}_{te_name}.fasta")

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



# above old







class Helpers:
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
        headers = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'expect value', 'bitscore']
        headers_dict = {i[0]: i[1] for i in enumerate(headers)}
        return pd.read_csv(filename, sep="\t", header=None).rename(headers_dict, axis=1)

    @staticmethod
    def save_df(df, filename):
        df.to_csv(filename, sep='\t', index=False, encoding='utf-8', header=False)







































import random
def rand_split(sequence, splits):
    for splitLen in splits:
        if splitLen > len(sequence):
            break
        yield sequence[:splitLen]
        sequence = sequence[splitLen::]

def rand_gen(int_list):
    while True:
        yield int(random.choice(int_list))

def cut_genome(genome_path, ont_path, coverage, filename):
    possibility_len = list(Helpers.get_len_fasta(ont_path).values())
    n = 0
    for i in range(coverage):
        print(i)
        for name, genome_sequence in Helpers.get_fasta(genome_path):
            print(name)
            with open(filename, 'a') as file:
                for seq in rand_split(genome_sequence, rand_gen(possibility_len)):
                    n = n + 1
                    file.write(f">{name}_{i}_{n}\n{seq}\n")









def mask_genome(te_name, te_path, genome_path, make_db=True):
    os.system(f"mkdir -p masked_genome")
    if make_db:
        os.system(f"makeblastdb -in {genome_path} -dbtype nucl -out genome/genome")
    os.system(f"blastn -num_threads 20 -outfmt 6 -query {te_path} -db genome/genome -out masked_genome/{te_name}_genome.bl -perc_identity 0.9")
    te_genome_bl = Helpers.read_bl(f"masked_genome/{te_name}_genome.bl")
    te_genome_bl = te_genome_bl[te_genome_bl.qseqid == te_name]
    te_genome_bl["start"] = te_genome_bl['sstart'].where(te_genome_bl['sstart'] <= te_genome_bl["send"], other=te_genome_bl['send'])
    te_genome_bl["stop"] = te_genome_bl['sstart'].where(te_genome_bl['sstart'] > te_genome_bl["send"], other=te_genome_bl['send'])
    # todo
    Helpers.save_df(te_genome_bl[["sseqid", "start", "stop"]], f"masked_genome/to_mask_genome_{te_name}.txt")
    plus_te_genome_bl = te_genome_bl[te_genome_bl.sstart < te_genome_bl.send]
    minus_te_genome_bl = te_genome_bl[te_genome_bl.sstart > te_genome_bl.send]

    merge_headers = ["sseqid", "sstart", "send"]
    Helpers.save_df(plus_te_genome_bl[merge_headers].sort_values(by=merge_headers), f"masked_genome/plus_to_mask_genome_{te_name}.bed")
    merge_headers = ["sseqid", "send", "sstart"]
    Helpers.save_df(minus_te_genome_bl[merge_headers].sort_values(by=merge_headers), f"masked_genome/minus_to_mask_genome_{te_name}.bed")

    os.system(f"bedtools merge -i masked_genome/plus_to_mask_genome_{te_name}.bed -d 200 > masked_genome/plus_to_mask_genome_{te_name}.bed.merged")
    os.system(f"bedtools merge -i masked_genome/minus_to_mask_genome_{te_name}.bed -d 200 > masked_genome/minus_to_mask_genome_{te_name}.bed.merged")
    minus_merged =  pd.read_csv(f"masked_genome/minus_to_mask_genome_{te_name}.bed.merged", sep="\t", header=None)
    plus_merged =  pd.read_csv(f"masked_genome/plus_to_mask_genome_{te_name}.bed.merged", sep="\t", header=None)
    merged_genome = pd.concat(
        [
            minus_merged[
                (minus_merged[2] - minus_merged[1] > from_te)
                & (minus_merged[2] - minus_merged[1] < to_te)
            ],
            plus_merged[
                (plus_merged[2] - plus_merged[1] > from_te)
                & (plus_merged[2] - plus_merged[1] < to_te)
            ],
        ]
    )
    merged_genome["length"] = merged_genome[2] - merged_genome[1]
    Helpers.save_df(merged_genome, f"masked_genome/to_mask_genome_{te_name}.bed.merged")
    os.system(f"bedtools maskfasta -fi {genome_path} -bed masked_genome/to_mask_genome_{te_name}.bed.merged -fo genome/{te_name}_masked_genome")
    return merged_genome


# first part
def create_subsequent_fasta(fasta, bed_df, file_base):
    # create subsequents with left and right flanks
    # site p or m
    num = 4000

    print(fasta)
    records = SeqIO.parse(fasta, 'fasta')
    record_dict = SeqIO.to_dict(records)
    for bed_rec in bed_df.itertuples():
        seq_id = bed_rec.sseqid
        record = record_dict[seq_id]

        te_start = bed_rec.sstart
        te_end = bed_rec.send

        left_start = te_start - num if te_start - num > 0 else 0
        left_end = te_start

        right_start = te_end
        right_end = te_end + num if te_end + num < len(record) else len(record)

        left_subseq_record = record[left_start:left_end]  # For the full record (with header)
        right_subseq_record = record[right_start:right_end]

        with open(f"first/left_{file_base}", 'a') as file:
            file.write(f">{record.id}\n{left_subseq_record.seq}\n")
        with open(f"first/right_{file_base}", 'a') as file:
            file.write(f">{record.id}\n{right_subseq_record.seq}\n")




"""
x(start_df, end_df, strain, ont_file) #cale 2 df
xx(start_record, end_df, strain) # jesli znajdzie dopasowuje prawo
xxx(start_record, end_record) # True False

"""


# second part


class SecondService:
    def __init__(self, filename, te_name):
        self.filename = filename
        self.te_name = te_name

        self.plus_te_start_records = {}
        self.plus_te_end_records = {}
        self.minus_te_start_records = {}
        self.minus_te_end_records = {}

    @property
    def used_ids(self):
        used_records = set()
        for _, records in self.plus_te_start_records.items():
            for record in records:
                used_records.add(record.qseqid)
        for _, records in self.minus_te_start_records.items():
            for record in records:
                used_records.add(record.qseqid)
        for _, records in self.plus_te_end_records.items():
            for record in records:
                used_records.add(record.qseqid)
        for _, records in self.minus_te_end_records.items():
            for record in records:
                used_records.add(record.qseqid)
        return used_records

    def add_start(self, record, ont_file, is_plus_strain):
        if is_plus_strain:
            if self.plus_te_start_records.get(ont_file) is None:
                self.plus_te_start_records[ont_file] = [record]
            else:
                self.plus_te_start_records[ont_file].append(record)
        else:
            if self.minus_te_start_records.get(ont_file) is None:
                self.minus_te_start_records[ont_file] = [record]
            else:
                self.minus_te_start_records[ont_file].append(record)

    def add_end(self, record, ont_file, is_plus_strain):
        if is_plus_strain:
            if self.plus_te_end_records.get(ont_file) is None:
                self.plus_te_end_records[ont_file] = [record]
            else:
                self.plus_te_end_records[ont_file].append(record)
        else:
            if self.minus_te_end_records.get(ont_file) is None:
                self.minus_te_end_records[ont_file] = [record]
            else:
                self.minus_te_end_records[ont_file].append(record)

    def get_start_field(self, is_plus_strain):
        return "send" if is_plus_strain else "sstart"

    def get_end_field(self, is_plus_strain):
        return "sstart" if is_plus_strain else "send"

    def if_match(self, start_record, end_record, is_plus_strain):
        start = getattr(start_record, self.get_start_field(is_plus_strain))
        end = getattr(end_record, self.get_end_field(is_plus_strain))
        return (end - start < to_te and end - start > from_te) or abs(end - start) < 200

    def match_rec_df(self, start_record, end_df, is_plus_strain):
        chromosome = start_record.sseqid
        r_df = end_df[(end_df.sseqid == chromosome)]
        for right_row in r_df.itertuples():
            if self.if_match(start_record, right_row, is_plus_strain):
                return right_row
        return None

    '''
    def match_not_ideal_pairs(self, left_df, right_df, is_plus_strain):
        start_df = left_df.sort_values(by=['sseqid', 'length', self.get_start_field(is_plus_strain)])
        end_df = right_df.sort_values(by=['sseqid', 'length', self.get_end_field(is_plus_strain)])
        if not len(start_df) or not len(end_df):
            return False
            if not len(start_df):
                if is_plus_strain:
                    self.plus_te_end_records.append(end_df.iloc[0])
                else:
                    self.minus_te_end_records.append(end_df.iloc[0])
                return True
            elif not len(end_df):
                if is_plus_strain:
                    self.plus_te_start_records.append(start_df.iloc[0])
                else:
                    self.minus_te_start_records.append(start_df.iloc[0])
                return True

            for chromosome in start_df.sseqid.unique():
                l_df =  start_df[(start_df.sseqid == chromosome)]
                r_df =  end_df[(end_df.sseqid == chromosome)]

        return False
    '''

    def match_ideal_pairs(self, left_df, right_df, ont_file, is_plus_strain):
        # add left and right if find
        # iterate chromosome
        # only one id
        start_df = left_df.sort_values(by=['sseqid', 'length', self.get_start_field(is_plus_strain)])
        end_df = right_df.sort_values(by=['sseqid', 'length', self.get_end_field(is_plus_strain)])
        if not len(start_df) or not len(end_df):
            return False
        for chromosome in start_df.sseqid.unique():
            l_df = start_df[(start_df.sseqid == chromosome)]
            r_df = end_df[(end_df.sseqid == chromosome)]
            for start_row in l_df.itertuples():
                matched_right_record = self.match_rec_df(start_row, r_df, is_plus_strain)
                if matched_right_record is not None:
                    if is_plus_strain:
                        self.add_start(start_row, ont_file, is_plus_strain)
                        self.add_end(matched_right_record, ont_file, is_plus_strain)
                    return True

        return False

    def save_to_file(self, record, is_start, is_plus_strain, ont_file):
        if record is not None:
            r = {
                "chromosome": record.sseqid,
                "strain": "plus" if is_plus_strain else "minus",
                "start": getattr(record, self.get_start_field(is_plus_strain)) if is_start else None,
                "stop": getattr(record, self.get_end_field(is_plus_strain)) if not is_start else None,
                "ont_id": f"{ont_file}:{record.qseqid}",
            }
            with open(self.filename, "a") as f_out:
                f_out.write(
                    f"{r['chromosome']}\t{r['strain']}\t{r['start']}\t{r['stop']}\t{r['ont_id']}\n"
                )

    def match_dfs(self, left_df, right_df, ont_file, is_plus_strain):
        # iterate ont_ids
        all_ids = set(left_df.qseqid.unique()).union(set(right_df.qseqid.unique()))
        # msstch 1
        for ont_id in all_ids:
            left = left_df[left_df.qseqid == ont_id]
            right = right_df[right_df.qseqid == ont_id]
            self.match_ideal_pairs(left, right, ont_file, is_plus_strain)

        # todo - now
        not_used_start_records = left_df[~left_df.qseqid.isin(self.used_ids)]
        not_used_end_records = right_df[~right_df.qseqid.isin(self.used_ids)]
        used_start_records = left_df[left_df.qseqid.isin(self.used_ids)]
        used_end_records = right_df[right_df.qseqid.isin(self.used_ids)]
        for not_used_start_record in not_used_start_records.itertuples():
            for used_start_record in used_start_records.itertuples():
                start1 = getattr(not_used_start_record, self.get_start_field(is_plus_strain))
                start2 = getattr(used_start_record, self.get_start_field(is_plus_strain))
                if abs(start1 - start2) < 300:
                    print('konik')

    # main function
    def xxx(self, left_blast, right_blast, ont_file):
        left_blast = left_blast[left_blast.length > 600]
        right_blast = right_blast[right_blast.length > 600]

        plus_left = left_blast[left_blast.sstart <= left_blast.send]
        plus_left["strain"] = "+"
        plus_right = right_blast[right_blast.sstart <= right_blast.send]
        plus_right["strain"] = "+"

        minus_left = left_blast[left_blast.sstart > left_blast.send]
        minus_left["strain"] = "-"
        minus_right = right_blast[right_blast.sstart > right_blast.send]
        minus_right["strain"] = "-"

        self.match_dfs(plus_left, plus_right, ont_file, True)
        self.match_dfs(minus_right, minus_left, ont_file, False)

    def ont_to_genome_blast(self, fasta_filename, blast_filename):
        return Blast.run(
            fasta_filename,
            f"genome/masked_genome",
            blast_filename,
        ).sort_values(['qseqid', 'bitscore'], ascending=False).groupby(['qseqid'], as_index=False).first()

    def blast(self, o, blast_db=True):
        left_m = self.ont_to_genome_blast(
            f"first/left_m_{o}_{self.te_name}.fasta",
            f"second/left_m_{o}_{self.te_name}_genome.bl",
        )
        right_m = self.ont_to_genome_blast(
            f"first/right_m_{o}_{self.te_name}.fasta",
            f"second/right_m_{o}_{self.te_name}_genome.bl",
        )
        left_p = self.ont_to_genome_blast(
            f"first/left_p_{o}_{self.te_name}.fasta",
            f"second/left_p_{o}_{self.te_name}_genome.bl",
        )
        right_p = self.ont_to_genome_blast(
            f"first/right_p_{o}_{self.te_name}.fasta",
            f"second/right_p_{o}_{self.te_name}_genome.bl",
        )
        return [
            [left_m, right_m], [left_p, right_p]
        ]

    def run(self, blast_db=True):

        raw_report_header = ["chromosome", "strain", "start", "stop", "ont_id"]
        header_to_save = '\t'.join(raw_report_header)
        with open(self.filename, "w") as f_out:
            f_out.write(
                f"{header_to_save}\n"
            )

        report = {}
        os.system("mkdir -p final")
        os.system("mkdir -p second")
        pd.options.mode.copy_on_write = True
        if blast_db:
            Blast.make_db(genome_path, "genome/masked_genome")
        for o in ont:
            for left_blast, right_blast in self.blast(o, blast_db):
                self.xxx(left_blast, right_blast, o)

        for ont_file, start_records in self.plus_te_start_records.items():
            for start_record in start_records:
                self.save_to_file(start_record, True, True, ont_file)
        for ont_file, end_records in self.plus_te_end_records.items():
            for end_record in end_records:
                self.save_to_file(end_record, False, True, ont_file)

        for ont_file, start_records in self.minus_te_start_records.items():
            for start_record in start_records:
                self.save_to_file(start_record, True, False, ont_file)
        for ont_file, end_records in self.minus_te_end_records.items():
            for end_record in end_records:
                self.save_to_file(end_record, False, False, ont_file)
        # report.to_csv(self.filename, sep='\t', index=False, encoding='utf-8', header=True)
        return pd.read_csv(self.filename, sep="\t")











def get_first_group(df, field):
    first_row = df.iloc[0]
    field_value = first_row[field]
    # group = df[(df[field] > field_value - 1700) & (df[field] < field_value + 1700)]
    # new_df = df[(df[field] <= field_value - 1700) | (df[field] >= field_value + 1700)]
    group = df[(df[field] > field_value - 1000) & (df[field] < field_value + 1000)]
    new_df = df[(df[field] <= field_value - 1000) | (df[field] >= field_value + 1000)]
    return new_df, group


def get_groups_dict(groups, field):
    group_dict = {}
    for group in groups:
        first_row = group.iloc[0]
        group_dict[first_row.chromosome, int(group[field].mode().median())] = group
    return group_dict


def get_groups(records, field):
    all_records = records[raw_report[field].notnull()]
    groups = []
    for chromosome in all_records.chromosome.unique():
        remaining = all_records[all_records.chromosome == chromosome]
        while len(remaining):
            g = get_first_group(remaining, field)
            remaining = g[0]
            groups.append(g[1])
    return groups


def if_match(start_chromosome, end_chromosome, start, end):
    if start_chromosome == end_chromosome:
        if (end - start > from_te and end - start < to_te) or (abs(end - start) < 200):
            return True
    return False


def create_record(chromosome, start, end, start_records, end_records):
    start_ids = set([ont_id for ont_id in start_records.ont_id])
    end_ids = set([ont_id for ont_id in end_records.ont_id])
    if start_ids.intersection(end_ids):
        pkt = 1
    elif set([x.split(":")[0] for x in start_ids]).intersection([x.split(":")[0] for x in end_ids]):
        pkt = 2
    elif start_ids and start_ids:
        pkt = 3
    else:
        pkt = 4
    return {
        "pkt": pkt,
        "chromosome": chromosome,
        "start": start,
        "stop": end,
        "length": end - start,
        "left_n": len(start_records),
        "left_rec": ",".join([ont_id for ont_id in start_records.ont_id]),
        "right_n": len(end_records),
        "right_rec": ",".join([ont_id for ont_id in end_records.ont_id]),
    }


def match(start_dict, end_dict):
    report_list = []
    start_remaining = start_dict.copy()
    end_remaining = end_dict.copy()
    for (start_chromosome, start), start_records in start_dict.items():
        end_dict = end_remaining.copy()
        for (end_chromosome, end), end_records in end_dict.items():
            if if_match(start_chromosome, end_chromosome, start, end):
                report_list.append(
                    create_record(
                        start_chromosome,
                        start,
                        end,
                        start_records,
                        end_records
                    )
                )
                del start_remaining[start_chromosome, start]
                del end_remaining[end_chromosome, end]
                continue
        # todo if not match?
    return report_list


def save_report(report_list, filename):
    report_header = ["pkt", "chromosome", "start", "stop", "length", "left_n", "left_rec", "right_n", "right_rec"]
    header_to_save = '\t'.join(report_header)
    with open(filename, "w") as f_out:
        f_out.write(
            f"{header_to_save}\n"
        )
    for r in report_list:
        with open(filename2, "a") as f_out:
            f_out.write(
                f"{r['pkt']}\t{r['chromosome']}\t{r['start']}\t{r['stop']}\t{r['length']}\t"
                f"{r['left_n']}\t{r['left_rec']}\t{r['right_n']}\t{r['right_rec']}\n"
            )


def save_first_point_report(filename):
    start_groups = get_groups(raw_report, "start")
    end_groups = get_groups(raw_report, "stop")

    start_dict = get_groups_dict(start_groups, "start")
    end_dict = get_groups_dict(end_groups, "stop")

    report_list = match(start_dict, end_dict)

    save_report(report_list, filename)

    report = pd.read_csv(filename, sep="\t").sort_values(by=['chromosome', 'pkt', 'start'])
    report.to_csv(filename, sep='\t', index=False, encoding='utf-8', header=True)
    return report


def _is_ref(df, chromosome, start, stop):
    return (
            (df.chromosome == chromosome)
            & (df.start < start + 500) & (df.start > start - 500)
            & (df.stop < stop + 500) & (df.stop > stop - 500)
    )


def is_ref(df, reference):
    is_reference = None
    for rec in masked_genome.itertuples():
        chromosome = rec[1]
        start = rec[2]
        stop = rec[3]
        part_is_reference = _is_ref(df, chromosome, start, stop)
        is_reference = is_reference | part_is_reference if is_reference is not None else part_is_reference
    return is_reference


def get_rec(x, o):
    records = []
    for r in x.split(","):
        rec = r.split(":")
        o_file = rec[0]
        o_id = rec[1]
        if o_file == o:
            records.append(o_id)
    return ",".join(records)


def get_num(x, o):
    records = []
    for r in x.split(","):
        rec = r.split(":")
        o_file = rec[0]
        o_id = rec[1]
        if o_file == o:
            records.append(o_id)
    return len(records)


def generate_second_report(filename_in, filename_out, reference):
    report_header = ["pkt", "chromosome", "start", "stop", "length", "is_ref"]
    df = pd.read_csv(filename2, sep="\t")
    for o in ont:
        df[f"left_{o}_num"] = df.left_rec.apply(
            get_num,
            args=[o]
        )
        df[f"left_{o}_rec"] = df.left_rec.apply(
            get_rec,
            args=[o]
        )
        report_header.append(f"left_{o}_num")
        report_header.append(f"left_{o}_rec")

        df[f"right_{o}_num"] = df.right_rec.apply(
            get_num,
            args=[o]
        )
        df[f"right_{o}_rec"] = df.right_rec.apply(
            get_rec,
            args=[o]
        )
        report_header.append(f"right_{o}_num")
        report_header.append(f"right_{o}_rec")

    is_reff = is_ref(df, masked_genome)
    df['is_ref'] = "nonref"
    df.loc[(df.length > 600), 'is_ref'] = "unknown"
    df.loc[is_reff, 'is_ref'] = "ref"

    print(report_header)
    report = df[report_header].sort_values(by=['chromosome', 'pkt', 'start'])
    report.to_csv(filename_out, sep='\t', index=False, encoding='utf-8', header=True)
    return report















from first_service import FirstService
from mask_genome import mask_genome
from second_service import SecondService
from second_service import ThirdService

if __name__ == "__main__":
    print("start")
    #okno = mask_genome(True)

    runn = FirstService()
    runn.run()

    """
    runn = SecondService()
    runn.run(False)

    runn = ThirdService()
    runn.run()
    """

    import ipdb; ipdb.set_trace()
    """
    masked_genome = mask_genome(te_name, te_path, genome_path, True)
    first_part(False, True, True)

    # cut_genome(genome_path, "ont/K10f.fasta", 20, 'cutted_genome.fasta')

    filename = "zizwa"
    filename2 = "moli2"
    filename3 = "wuxi3"

    service = SecondService(filename, te_name)
    raw_report = service.run(False)

    first_report = save_first_point_report(filename2)
    xxx = generate_second_report(filename2, filename3, masked_genome)
    """

    #old
    """
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
    """

