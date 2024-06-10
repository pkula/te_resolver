import logging
import os
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from config import Config
from blast import Blast
from helpers import Helpers


class FirstService:
    def __init__(self):
        self.config = Config()

    #checked

    def _first_part(self, ont_base, te_ont_bl, is_plus_strain):
        o = self.config.get_te_ont_bl_path(ont_base)
        x = "p" if is_plus_strain else "m"

        basename = f"{x}_{o.stem}"
        bed_filename = self.config.first_bed_path / f'{basename}.bed'
        merge_filename = self.config.first_bed_path / f"{basename}.bed.merged"
        txt_filename = self.config.first_bed_path / f'{basename}_seq_list.txt'
        raw_subseq_filename = self.config.first_subseq_path / f"raw_{basename}_subseq.fasta"
        subseq_filename = self.config.first_subseq_path / f"{basename}_subseq.fasta"

        fasta_file = self.config.get_ont_filepath_from_ont_base(ont_base)

        start_field = "sstart" if is_plus_strain else "send"
        end_field = "send" if is_plus_strain else "sstart"

        te_ont_bed = te_ont_bl[['sseqid', start_field, end_field, 'qseqid']].sort_values(
            by=['sseqid', start_field, end_field])
        Helpers.save_df(te_ont_bed, bed_filename)
        os.system(
            f"bedtools merge -i {bed_filename} -d 1500 -c 4 -o collapse  > {merge_filename}")
        merged = pd.read_csv(
            f"{merge_filename}", sep="\t", header=None, names=['sseqid', 'sstart', 'send', 'qseqid']
        )
        filtered_bed = merged[(merged.send - merged.sstart) > self.config.from_te]
        Helpers.save_df(pd.DataFrame({'id': filtered_bed['sseqid'].unique()}), txt_filename)
        os.system(
            f"seqtk subseq {fasta_file} {txt_filename} > {raw_subseq_filename}")
        self.create_subsequent_fasta(raw_subseq_filename, filtered_bed, subseq_filename)


    def first_part(self):
        for ont_base in self.config.ont_bases:
            te_ont_bl = Helpers.read_bl(self.config.get_te_ont_bl_path(ont_base))
            te_ont_bl = te_ont_bl[te_ont_bl.qseqid.str.contains(self.config.te_name)]

            self._first_part(ont_base, te_ont_bl[te_ont_bl.sstart <= te_ont_bl.send], True)
            self._first_part(ont_base, te_ont_bl[te_ont_bl.sstart > te_ont_bl.send], False)


    def create_subsequent_fasta(self, fasta, bed_df, filename):
        # create subsequents with left and right flanks
        # site p or m

        # todo filenames
        parent = filename.parent
        f = filename.name
        with open(parent / f"{'left_' + f}", 'w'):
            pass
        with open(parent / f"{'right_' + f}", 'w'):
            pass

        records = SeqIO.parse(fasta, 'fasta')
        record_dict = SeqIO.to_dict(records)
        for bed_rec in bed_df.itertuples():
            seq_id = bed_rec.sseqid
            record = record_dict[seq_id]

            te_start = bed_rec.sstart
            te_end = bed_rec.send

            left_start = te_start - self.config.flanks_len if te_start - self.config.flanks_len > 0 else 0
            left_end = te_start

            right_start = te_end
            right_end = te_end + self.config.flanks_len if te_end + self.config.flanks_len < len(record) else len(record)

            left_subseq_record = record[left_start:left_end]  # For the full record (with header)
            right_subseq_record = record[right_start:right_end]

            with open(parent / f"{'left_' + f}", 'a') as file:
                file.write(f">{record.id}\n{left_subseq_record.seq}\n")
            with open(parent / f"{'right_' + f}", 'a') as file:
                file.write(f">{record.id}\n{right_subseq_record.seq}\n")

    def convert_fq(self):
        logging.info("start converting pass from fastaq to fasta")
        for f in self.config.ont_path.glob("*.fastq.gz"):
            self.from_fq_to_fa(f)
        logging.info("end converting pass to fasta")

    @staticmethod
    def from_fq_to_fa(fq):
        # change fq.gz into fasta (for blast)
        fq = Path(fq).absolute()
        basename = Path(fq).name.split(".")[0]
        os.system(f"seqtk seq -a {fq} > {fq.parent / basename}.fasta")

    def run(self, fq=True, create_db=True, do_blast=True):
        print("start first service")
        if fq:
            print("convert fq")
            self.convert_fq()
        if create_db:
            print("make blast")
            Blast.make_ont_db(self.config)
        if do_blast:
            print("do blast")
            Blast.run_te_ont(self.config)
        self.first_part()
        print("end first service")
