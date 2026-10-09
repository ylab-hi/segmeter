import os
import subprocess
import platform
import re

def sort_BED(infile, outfile):
    """Sorted copy of a BED file (unmeasured); `LC_ALL=C` so that the order (chromosome names, scaffolds) is the same on every
    machine as in the containers, which a seeded simulation needs (#15)"""
    with open(outfile, 'w') as out:
        subprocess.run(["sort", "-k1,1", "-k2,2n", "-k3,3n", str(infile)], stdout=out, env={**os.environ, "LC_ALL": "C"})

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

def get_time_rss_per_mb():
    """`/usr/bin/time` reports the RSS in bytes on macOS (`-l`) and in kilobytes on Linux (GNU time `-v`)"""
    if get_os() == "macos":
        return 1024 * 1024
    else: # linux
        return 1024

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


def chrom_lengths(paths):
    """Largest end per chromosome over the BED files, skipping blank, `#`, `track` and `browser` lines."""
    lengths = {}
    for path in paths:
        with open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                try:
                    chrom, _, end = line.split("\t")[:3]
                    lengths[chrom] = max(lengths.get(chrom, 0), int(end))
                except ValueError:
                    raise ValueError(f"{path}: not a tab-separated BED line: {line!r}") from None
    return lengths


# basic query types and their group: interval (overlaps the reference) or gap (does not)
BASIC_QUERIES = {
    "perfect": "interval", "5p-partial": "interval", "3p-partial": "interval", "enclosed": "interval", "contained": "interval",
    "perfect-gap": "gap", "left-adjacent-gap": "gap", "right-adjacent-gap": "gap", "mid-gap1": "gap", "mid-gap2": "gap",
}

def get_query_group(query):
    """returns the group of a basic query type: interval (overlaps the reference) or gap (does not)"""
    if query not in BASIC_QUERIES:
        raise ValueError("Query type not supported")
    return BASIC_QUERIES[query]
