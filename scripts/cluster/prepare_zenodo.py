#!/usr/bin/env python3
"""Download, verify and extract the fixed published simulation archive."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import tarfile
import urllib.request

URL = 'https://zenodo.org/api/records/14880992/files/simdata.tar.gz/content'
CHECKSUM = '75b3c767def3813aa790d01e9ae9239a'


def verify(path):
    digest = hashlib.md5()
    with path.open('rb') as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b''):
            digest.update(block)
    if digest.hexdigest() != CHECKSUM:
        raise ValueError(f'Zenodo checksum mismatch: {path}; remove the bad archive and retry')


def extract(archive, directory):
    prefix = PurePosixPath('simdata/sim/sim_001/BED')
    with tarfile.open(archive, 'r:gz') as tar:
        for member in tar:
            path = PurePosixPath(member.name)
            if path.is_absolute() or '..' in path.parts or member.issym() or member.islnk():
                raise ValueError(f'Unsafe archive entry: {member.name}')
            if path != prefix and prefix not in path.parents:
                continue
            if member.isdir():
                continue
            if not member.isfile():
                raise ValueError(f'Unsupported archive entry: {member.name}')
            if path.name.startswith('.'):
                continue
            target = directory.joinpath(*path.parts)
            target.parent.mkdir(parents=True, exist_ok=True)
            with tar.extractfile(member) as src, target.open('wb') as dst:
                shutil.copyfileobj(src, dst)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, required=True)
    args = parser.parse_args()
    directory = args.directory.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    archive = directory / 'simdata.tar.gz'
    if not archive.exists():
        partial = directory / 'simdata.tar.gz.part'
        print('Downloading Zenodo simulation archive (216 MB)...', flush=True)
        with urllib.request.urlopen(URL, timeout=120) as src, partial.open('wb') as dst:
            shutil.copyfileobj(src, dst)
        verify(partial)
        partial.rename(archive)
    verify(archive)
    marker = directory / 'simdata/zenodo.json'
    if not marker.exists():
        print('Extracting verified simulation data...', flush=True)
        extract(archive, directory)
        marker.write_text(json.dumps({'record': '14880992', 'url': URL, 'md5': CHECKSUM}, indent=2))
    print(f'Zenodo data ready: {directory / "simdata"}')


if __name__ == '__main__':
    main()
