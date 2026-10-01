import subprocess
import platform
import re

def sort_BED(infile, outfile):
    with open(outfile, 'w') as out:
        subprocess.run(["sort", "-k1,1", "-k2,2n", "-k3,3n", str(infile)], stdout=out)

def file_linecounter(filepath):
    with open(filepath, "rb") as file:
        return sum(1 for _ in file)

def get_os():
    """returns the operating system"""
    if platform.system() == "Darwin":
        return "macos"
    elif platform.system() == "Linux":
        return "linux"
    else:
        raise ValueError("Operating system not supported")

def get_time_rss_label():
    if get_os() == "macos":
        return "maximum resident set size"
    else: # linux
        return "Maximum resident set size (kbytes)"

def get_time_verbose_flag():
    if get_os() == "macos":
        return "-l"
    else:
        return "-v"

def get_rss_from_stderr(stderr_output, rss_label):
    for line in stderr_output.split("\n"):
        if rss_label in line:
            # Extract the numerical value from the line
            match = re.search(r"(\d+)", line)
            if match:
                return int(match.group(1))
    return -1


def chrom_sort_key(name):
    """Chromosome order of genomap 0.2.6 (`chromosome_probe`), the map type of granges: numbers first,
    then X, Y, M, Z, W, O, then the rest. granges 0.2.2 labels the query trees in this order while it
    reads the genome file in file order, so a genome file in any other order makes it compare the
    wrong chromosomes."""
    name = name[3:] if name.startswith("chr") else name
    letters = {"X": 2, "Y": 3, "M": 4, "MT": 4, "Mt": 4, "Z": 5, "W": 6, "O": 7}
    arm = {"L": 1, "R": 2}.get(name[-1:], 0)
    number = name[:-1] if arm else name
    if number.isdigit():
        return (1, int(number), "", arm)
    return (letters.get(name, 8), 0, name, 0)


def get_query_group(datatype, query):
    """returns the query group (e.g., interval or gap) based on the datatype and query"""

    if datatype == "basic":
        if query in ["perfect", "5p-partial", "3p-partial", "enclosed", "contained"]:
            return "interval"
        elif query in ["perfect-gap", "left-adjacent-gap", "right-adjacent-gap", "mid-gap1", "mid-gap2"]:
            return "gap"
        else:
            raise ValueError("Query type not supported")
    elif datatype == "complex":
        return "undefined" # todo
