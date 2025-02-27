import subprocess
import platform
import re

def sort_BED(infile, outfile):
    with open(outfile, 'w') as out:
        subprocess.run(["sort", "-k1,1", "-k2,2n", "-k3,3n", str(infile)], stdout=out)

def file_linecounter(filepath):
    total = 0
    with open(filepath, "rb") as file:
        total += 1
    return total

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


def save_index_time(self, intvlnum, index_time, filename):
    fh = open(filename, "w")
    fh.write("intvlnum\ttime(s)\n")
    for key, value in index_time.items():
        fh.write(f"{key}\t{value}\n")
    fh.close()


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
