from pathlib import Path
import logging
from sqlite3.dbapi2 import Date
import subprocess

def sort_BED(infile, outfile):
    with open(outfile, 'w') as out:
        subprocess.run(["sort", "-k1,1", "-k2,2n", "-k3,3n", str(infile)], stdout=out)

def file_linecounter(filepath):
    try:
        with open(filepath, "rb") as file:
            return(sum(1 for line in file))
    except FileNotFoundError:
        print(f"File {filepath} not found.")
    except Exception as e:
        print(f"Error: {e}")


def get_log_level(level):
    loglevel = {
        "DEBUG": logging.DEBUG,
        "INFO": logging.INFO,
        "WARNING": logging.WARNING,
        "ERROR": logging.ERROR,
        "CRITICAL": logging.CRITICAL
    }
    return loglevel[level]

def setup_logger(logfile, loglevel, tool):
    logger = logging.getLogger(f"{tool}")
    logger.setLevel(get_log_level(loglevel))

    if not logger.hasHandlers():
        fh = logging.FileHandler(logfile)
        fh.setLevel(get_log_level(loglevel))
        fh.setFormatter(logging.Formatter(
            f"%(asctime)s:%(name)s - %(message)s (%(levelname)s)", datefmt='%d-%m-%y %H:%M:%S'))
        logger.addHandler(fh)

        # Console handler
        sh = logging.StreamHandler()
        sh.setLevel(get_log_level(loglevel))
        sh.setFormatter(logging.Formatter(
            f"%(asctime)s:%(name)s - %(message)s (%(levelname)s)", datefmt='%d-%m-%y %H:%M:%S'))
        logger.addHandler(sh)

    return logger
