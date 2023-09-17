import os
import argparse



if __name__ == "__main__":
    # prepare parser
    parser = argparse.ArgumentParser(description="TE pipeline resolver")
    parser.add_argument('-w', action='store_true')
    args = parser.parse_args()
