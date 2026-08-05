import csv
import gzip
import hashlib
import io
import json
import shlex
import subprocess
from pathlib import Path


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def gzip_bytes(data):
    bio = io.BytesIO()
    with gzip.GzipFile(fileobj=bio, mode="wb", filename="", mtime=0, compresslevel=6) as gz:
        gz.write(data)
    return bio.getvalue()


def write_gzip(path, text):
    Path(path).write_bytes(gzip_bytes(text.encode()))


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")


def write_tsv(path, fields, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(rows)


def printable(argv):
    return shlex.join([str(x) for x in argv])


def run(argv, commands=None, stdout=None):
    if commands is not None:
        commands.append(printable(argv))
    subprocess.run([str(x) for x in argv], check=True, stdout=stdout)
