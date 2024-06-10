import argparse
import logging
import os
from pathlib import Path
import pandas as pd
from bioinfokit.analys import Fasta
from Bio import SeqIO
import numpy as np
from src.blast import Blast
from src.helpers import Helpers
from src.config import Config
from src.mask_genome import mask_genome


class SecondService:
    def __init__(self):
        self.config = Config()
        self.te_name = self.config.te_name
        self.masked_database = self.config.masked_genome_path / f"masked_{Helpers.get_filebase(self.config.genome_filepath)}"

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
        return (end - start < self.config.to_te and end - start > self.config.from_te) or abs(end - start) < self.config.max_nonref_len

    def match_rec_df(self, start_record, end_df, is_plus_strain):
        chromosome = start_record.sseqid
        r_df = end_df[(end_df.sseqid == chromosome)]
        for right_row in r_df.itertuples():
            if self.if_match(start_record, right_row, is_plus_strain):
                return right_row
        return None

    def match_pairs(self, left_df, right_df, ont_file, is_plus_strain):
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
            with open(self.config.filtered_records_filepath, "a") as f_out:
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
            self.match_pairs(left, right, ont_file, is_plus_strain)
        # msstch 2
        is_the_same_n = False
        while not is_the_same_n:
            n = len(self.used_ids)
            not_used_start_records = left_df[~left_df.qseqid.isin(self.used_ids)]
            not_used_end_records = right_df[~right_df.qseqid.isin(self.used_ids)]
            self.match_pairs(not_used_start_records, not_used_end_records, ont_file, is_plus_strain)
            if len(self.used_ids) == n:
                is_the_same_n = True

    # main function
    def add_start_and_end_te(self, left_blast, right_blast, ont_file):
        left_blast = left_blast[left_blast.length > self.config.group_len]
        right_blast = right_blast[right_blast.length > self.config.group_len]

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
            self.masked_database,
            blast_filename,
            self.config.blast_threads
        ).sort_values(['qseqid', 'bitscore'], ascending=False).groupby(['qseqid'], as_index=False).first()

    def blast(self, o, blast_db=True):
        left_m = self.ont_to_genome_blast(
            self.config.first_subseq_path / f"left_m_TE_{self.config.genome_filepath.stem}_{o}_{self.te_name}_subseq.fasta",
            self.config.second_path / f"left_m_{self.config.genome_filepath.stem}_{o}_{self.te_name}_genome.bl",
        )
        right_m = self.ont_to_genome_blast(
            self.config.first_subseq_path / f"right_m_TE_{self.config.genome_filepath.stem}_{o}_{self.te_name}_subseq.fasta",
            self.config.second_path / f"right_m_{self.config.genome_filepath.stem}_{o}_{self.te_name}_genome.bl",
        )
        left_p = self.ont_to_genome_blast(
            self.config.first_subseq_path / f"left_p_TE_{self.config.genome_filepath.stem}_{o}_{self.te_name}_subseq.fasta",
            self.config.second_path / f"left_p_{self.config.genome_filepath.stem}_{o}_{self.te_name}_genome.bl",
        )
        right_p = self.ont_to_genome_blast(
            self.config.first_subseq_path / f"right_p_TE_{self.config.genome_filepath.stem}_{o}_{self.te_name}_subseq.fasta",
            self.config.second_path / f"right_p_{self.config.genome_filepath.stem}_{o}_{self.te_name}_genome.bl",
        )
        return [
            [left_m, right_m], [left_p, right_p]
        ]

    def run(self, blast_db=True):
        # return raw report

        raw_report_header = ["chromosome", "strain", "start", "stop", "ont_id"]
        header_to_save = '\t'.join(raw_report_header)
        with open(self.config.filtered_records_filepath, "w") as f_out:
            f_out.write(
                f"{header_to_save}\n"
            )

        pd.options.mode.copy_on_write = True
        if blast_db:
            Blast.make_db(self.config.genome_filepath, self.masked_database)
        for o in self.config.ont_bases:
            for left_blast, right_blast in self.blast(o, blast_db):
                self.add_start_and_end_te(left_blast, right_blast, o)

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
        return pd.read_csv(self.config.filtered_records_filepath, sep="\t")



class ThirdService:
    def __init__(self):
        self.config = Config()

    def get_first_group(self, df, field):
        first_row = df.iloc[0]
        field_value = first_row[field]
        group = df[(df[field] > field_value - self.config.group_len) & (df[field] < field_value + self.config.group_len)]
        new_df = df[(df[field] <= field_value - self.config.group_len) | (df[field] >= field_value + self.config.group_len)]
        return new_df, group


    def get_groups_dict(self, groups, field):
        group_dict = {}
        for group in groups:
            first_row = group.iloc[0]
            group_dict[first_row.chromosome, int(group[field].mode().median())] = group
        return group_dict


    def get_groups(self, records, field):
        all_records = records[self.raw_report[field].notnull()]
        groups = []
        for chromosome in all_records.chromosome.unique():
            remaining = all_records[all_records.chromosome == chromosome]
            while len(remaining):
                g = self.get_first_group(remaining, field)
                remaining = g[0]
                groups.append(g[1])
        return groups


    def if_match(self, start_chromosome, end_chromosome, start, end):
        if start_chromosome == end_chromosome:
            if (end - start > self.config.from_te and end - start < self.config.to_te) or (abs(end - start) < self.config.max_nonref_len):
                return True
        return False


    def create_record(self, chromosome, start, end, start_records, end_records):
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


    def match_end_dict(self, start_chromosome, start, start_records, end_dict):
        for (end_chromosome, end), end_records in end_dict.items():
            if self.if_match(start_chromosome, end_chromosome, start, end):
                return (end_chromosome, end),  self.create_record(
                    start_chromosome,
                    start,
                    end,
                    start_records,
                    end_records
                )

    def match_dicts(self, start_dict, end_dict):
        report_list = []
        start_remaining = start_dict.copy()
        end_remaining = end_dict.copy()
        for (start_chromosome, start), start_records in start_dict.items():
            end_dict = end_remaining.copy()
            record = self.match_end_dict(start_chromosome, start, start_records, end_dict)
            if record:
                report_list.append(
                    record[1]
                )
                del start_remaining[start_chromosome, start]
                del end_remaining[record[0][0], record[0][1]]
            # todo if not match?
        return report_list


    def save_report(self, report_list):
        report_header = ["pkt", "chromosome", "start", "stop", "length", "left_n", "left_rec", "right_n", "right_rec"]
        header_to_save = '\t'.join(report_header)
        with open(self.config.raw_report_filepath, "w") as f_out:
            f_out.write(
                f"{header_to_save}\n"
            )
        for r in report_list:
            with open(self.config.raw_report_filepath, "a") as f_out:
                f_out.write(
                    f"{r['pkt']}\t{r['chromosome']}\t{r['start']}\t{r['stop']}\t{r['length']}\t"
                    f"{r['left_n']}\t{r['left_rec']}\t{r['right_n']}\t{r['right_rec']}\n"
                )


    def save_first_point_report(self):

        start_groups = self.get_groups(self.raw_report, "start")
        end_groups = self.get_groups(self.raw_report, "stop")

        start_dict = self.get_groups_dict(start_groups, "start")
        end_dict = self.get_groups_dict(end_groups, "stop")

        report_list = self.match_dicts(start_dict, end_dict)

        self.save_report(report_list)

        report = pd.read_csv(self.config.raw_report_filepath, sep="\t").sort_values(by=['chromosome', 'pkt', 'start'])
        report.to_csv(self.config.raw_report_filepath, sep='\t', index=False, encoding='utf-8', header=True)
        return report


    def _is_ref(self, df, chromosome, start, stop):
        return (
                (df.chromosome == chromosome)
                & (df.start < start + self.config.group_len) & (df.start > start - self.config.group_len)
                & (df.stop < stop + self.config.group_len) & (df.stop > stop - self.config.group_len)
        )


    def is_ref(self, df, reference):
        is_reference = None
        for rec in reference.itertuples():
            chromosome = rec[1]
            start = rec[2]
            stop = rec[3]
            part_is_reference = self._is_ref(df, chromosome, start, stop)
            is_reference = is_reference | part_is_reference if is_reference is not None else part_is_reference
        return is_reference


    def get_rec(self, x, o):
        records = []
        for r in x.split(","):
            rec = r.split(":")
            o_file = rec[0]
            o_id = rec[1]
            if o_file == o:
                records.append(o_id)
        return ",".join(records)


    def get_num(self, x, o):
        records = []
        for r in x.split(","):
            rec = r.split(":")
            o_file = rec[0]
            o_id = rec[1]
            if o_file == o:
                records.append(o_id)
        return len(records)


    def generate_second_report(self, reference):
        report_header = ["pkt", "chromosome", "start", "stop", "length", "is_ref"]
        df = pd.read_csv(self.config.raw_report_filepath, sep="\t")
        for o in self.config.ont_bases:
            df[f"left_{o}_num"] = df.left_rec.apply(
                self.get_num,
                args=[o]
            )
            df[f"left_{o}_rec"] = df.left_rec.apply(
                self.get_rec,
                args=[o]
            )
            report_header.append(f"left_{o}_num")
            report_header.append(f"left_{o}_rec")

            df[f"right_{o}_num"] = df.right_rec.apply(
                self.get_num,
                args=[o]
            )
            df[f"right_{o}_rec"] = df.right_rec.apply(
                self.get_rec,
                args=[o]
            )
            report_header.append(f"right_{o}_num")
            report_header.append(f"right_{o}_rec")

        is_reff = self.is_ref(df, reference)
        df['is_ref'] = "nonref"
        df.loc[(df.length > self.config.group_len), 'is_ref'] = "unknown"
        df.loc[is_reff, 'is_ref'] = "ref"

        report = df[report_header].sort_values(by=['chromosome', 'pkt', 'start'])
        report.to_csv(self.config.final_report_filepath, sep='\t', index=False, encoding='utf-8', header=True)
        return report

    def run(self):
        self.raw_report = pd.read_csv(self.config.filtered_records_filepath, sep="\t")
        reference = mask_genome()
        first_report = self.save_first_point_report()
        second_report = self.generate_second_report(reference)
        return first_report, second_report

