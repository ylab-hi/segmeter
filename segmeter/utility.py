from pathlib import Path
import subprocess

def sort_BED(infile, outfile):
    with open(outfile, 'w') as out:
        subprocess.run(["sort", "-k1,1", "-k2,2n", "-k3,3n", str(infile)], stdout=out)

def file_linecounter(filepath):
    total = 0
    with open(filepath, "rb") as file:
        total += 1
    return total





    # try:
    #     with open(filepath, "rb") as file:
    #         return(sum(1 for line in file))
    # except FileNotFoundError:
    #     print(f"File {filepath} not found.")
    # except Exception as e:
    #     print(f"Error: {e}")
