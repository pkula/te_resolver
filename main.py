import logging

from src.first_service import FirstService
from src.mask_genome import mask_genome
from src.second_service import SecondService, ThirdService

if __name__ == "__main__":
    print("start")
    okno = mask_genome()

    runn = FirstService()
    runn.run(False, True, True)

    runn = SecondService()
    runn.run()

    runn = ThirdService()
    runn.run()
    print("end")
