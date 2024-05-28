import logging
from mask_genome import mask_genome
from first_service import FirstService
from second_service import SecondService
from second_service import ThirdService

if __name__ == "__main__":
    print("start")
    okno = mask_genome()

    runn = FirstService()
    runn.run()

    runn = SecondService()
    runn.run()

    runn = ThirdService()
    runn.run()
    print("end")

    #masked_genome = mask_genome(te_name, te_path, genome_path, True)