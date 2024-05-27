import logging
import os
from pathlib import Path
import pandas as pd
from helpers import Helpers
from config import Config


def mask_genome(make_db=True):
    config = Config()
    te_name = config.te_name

    genome_base = Helpers.get_filebase(config.genome_filepath)
    filebase = f"{te_name}_{genome_base}"
    database_filepath = config.masked_genome_path / Helpers.get_filebase(config.genome_filepath)

    os.system(f"mkdir -p masked_genome")
    if make_db:
        os.system(f"makeblastdb -in {config.genome_filepath} -dbtype nucl -out {database_filepath}")
    os.system(f"blastn -num_threads 20 -outfmt 6 -query {config.te_filepath} -db {database_filepath} -out masked_genome/{filebase}.bl -perc_identity 0.9")
    te_genome_bl = Helpers.read_bl(f"masked_genome/{filebase}.bl")
    te_genome_bl = te_genome_bl[te_genome_bl.qseqid == te_name]
    te_genome_bl["start"] = te_genome_bl['sstart'].where(te_genome_bl['sstart'] <= te_genome_bl["send"], other=te_genome_bl['send'])
    te_genome_bl["stop"] = te_genome_bl['sstart'].where(te_genome_bl['sstart'] > te_genome_bl["send"], other=te_genome_bl['send'])
    # todo
    Helpers.save_df(te_genome_bl[["sseqid", "start", "stop"]], f"masked_genome/{filebase}.txt")
    plus_te_genome_bl = te_genome_bl[te_genome_bl.sstart < te_genome_bl.send]
    minus_te_genome_bl = te_genome_bl[te_genome_bl.sstart > te_genome_bl.send]

    merge_headers = ["sseqid", "sstart", "send"]
    Helpers.save_df(plus_te_genome_bl[merge_headers].sort_values(by=merge_headers), f"masked_genome/plus_{filebase}.bed")
    merge_headers = ["sseqid", "send", "sstart"]
    Helpers.save_df(minus_te_genome_bl[merge_headers].sort_values(by=merge_headers), f"masked_genome/minus_{filebase}.bed")

    os.system(f"bedtools merge -i masked_genome/plus_{filebase}.bed -d 200 > masked_genome/plus_{filebase}.bed.merged")
    os.system(f"bedtools merge -i masked_genome/minus_{filebase}.bed -d 200 > masked_genome/minus_{filebase}.bed.merged")

    plus_merged =  pd.read_csv(f"masked_genome/plus_{filebase}.bed.merged", sep="\t", header=None)
    minus_merged =  pd.read_csv(f"masked_genome/minus_{filebase}.bed.merged", sep="\t", header=None)
    merged_genome = pd.concat(
        [
            minus_merged[
                (minus_merged[2] - minus_merged[1] > config.from_te)
                & (minus_merged[2] - minus_merged[1] < config.to_te)
            ],
            plus_merged[
                (plus_merged[2] - plus_merged[1] > config.from_te)
                & (plus_merged[2] - plus_merged[1] < config.to_te)
            ],
        ]
    )
    merged_genome["length"] = merged_genome[2] - merged_genome[1]
    Helpers.save_df(merged_genome, f"masked_genome/{filebase}.bed.merged")
    os.system(f"bedtools maskfasta -fi {config.genome_filepath} -bed masked_genome/{filebase}.bed.merged -fo genome/{filebase}_masked.fasta")
    return merged_genome

