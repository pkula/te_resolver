import logging

from src.first_service import FirstService
from src.mask_genome import mask_genome
from src.second_service import SecondService, ThirdService
from src.config import Config
from src.cutter import cut_genome

if __name__ == "__main__":
    print("start")
    config = Config()
    if config.script == "finder":
        runn = FirstService()
        runn.run(config.convert_fq, config.create_db, config.do_blast)

        runn = SecondService()
        runn.run()

        runn = ThirdService()
        runn.run()
    elif config.script == "cutter":
        cut_genome(config.genome_filepath, config.cutter_sequencing_filename, config.cutter_coverage if config.cutter_coverage else 5, config.cutter_out_filename, config.cutter_from_n if config.cutter_from_n else 0)
    elif config.script == "masker":
        mask_genome()
    print("end")
