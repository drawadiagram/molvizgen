#!/usr/bin/env python3
"""
Shared I/O for the "candidate manifest" contract passed between FIND, FILTER,
and (optionally) GENERATE steps of a SELECT|GENERATE|ASSEMBLE pipeline.

A manifest is a JSON object:
    {"candidates": [
        {"id": "4e34", "pdb_path": "/abs/path/4e34.pdb",
         "chain_domain": "A", "chain_peptide": "B",
         ... any other scoring/metadata fields ...},
        ...
    ]}

FIND steps produce manifests (with or without score fields attached).
FILTER steps consume a manifest and produce a smaller (or reordered) manifest.
Every record must carry an absolute `pdb_path` so downstream steps never need
to guess file locations from an id + a directory convention.
"""
import json
import os


def write_manifest(path, candidates):
    """Write a list of candidate dicts to `path` as a manifest JSON object."""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w") as f:
        json.dump({"candidates": candidates}, f, indent=2)


def read_manifest(path):
    """Read a manifest JSON object (or file) and return its candidate list."""
    with open(path) as f:
        data = json.load(f)
    return data["candidates"]


def read_manifest_stdin_or_path(path_or_dash):
    """Read a manifest from a file path, or from stdin if given '-'."""
    if path_or_dash == "-":
        import sys
        data = json.load(sys.stdin)
        return data["candidates"]
    return read_manifest(path_or_dash)
