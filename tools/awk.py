import argparse

def main():
    options = parse_arguments()

def parse_arguments():
    parser = argparse.ArgumentParser(description="Interval overlap using awk")
    parser.add_argument("-q", "--query", type=str, required=True, help="Query intervals BED file")
    parser.add_argument("-t", "--target", type=str, required=True, help="Reference intervals BED file")


main()
