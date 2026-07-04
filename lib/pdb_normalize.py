#!/usr/bin/env python3
"""
Normalize non-standard, multi-character PDB chain identifiers into valid,
standards-compliant single-character chain IDs.

Why this exists: the strict PDB fixed-column format reserves exactly one
column (22) for the chain ID. Some pipelines (e.g. Boltz) write
multi-character tokens there (observed: "pdz"/"pep") using free-form,
whitespace-separated fields instead of fixed columns. Readers that parse
strict fixed columns (PyMOL included) then silently truncate the chain ID to
its first character and misread every field to its right (residue number,
coordinates, ...) — for "pdz"/"pep" both truncate to the same letter "p",
making chain-based selection impossible and quietly wrong rather than
erroring.

This module re-tokenizes ATOM/HETATM/TER lines by whitespace (the only
reliable way to recover the intended fields from such a file) and re-emits
them in correct fixed-column PDB format with a caller-supplied single-letter
chain remapping.
"""
import os


def _format_atom_name(name):
    if len(name) >= 4:
        return name[:4]
    return " " + name.ljust(3)


def _format_atom_line(tokens, record, chain_map):
    # record, serial, name, resName, chain, resSeq, x, y, z, occ, temp, element
    _record, serial, name, resname, chain, resseq, x, y, z, occ, temp, element = tokens
    new_chain = chain_map[chain]
    return (
        f"{record:<6}"
        f"{int(serial):>5} "
        f"{_format_atom_name(name)}"
        f" "
        f"{resname:>3} "
        f"{new_chain:1}"
        f"{int(resseq):>4}"
        f"    "
        f"{float(x):>8.3f}{float(y):>8.3f}{float(z):>8.3f}"
        f"{float(occ):>6.2f}{float(temp):>6.2f}"
        f"          {element:>2}\n"
    )


def _format_ter_line(tokens, chain_map):
    _record, serial, resname, chain, resseq = tokens
    new_chain = chain_map[chain]
    return (
        f"{'TER':<6}"
        f"{int(serial):>5} "
        f"     "
        f"{resname:>3} "
        f"{new_chain:1}"
        f"{int(resseq):>4}\n"
    )


def detect_chain_tokens(pdb_path):
    """Return the distinct whitespace-tokenized chain ids in a PDB file's
    ATOM records, in first-seen order."""
    seen = []
    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            chain = parts[4]
            if chain not in seen:
                seen.append(chain)
    return seen


def normalize_pdb_chains(in_path, out_path, chain_map):
    """Rewrite `in_path` to `out_path`, remapping whitespace-tokenized chain
    ids per `chain_map` (e.g. {"pdz": "A", "pep": "B"}) into standard,
    single-character, fixed-column PDB chain IDs. Lines whose record type
    isn't ATOM/HETATM/TER are copied through unchanged. Raises KeyError if a
    chain token is encountered that isn't in `chain_map`.
    """
    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    with open(in_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                parts = line.split()
                if len(parts) < 12:
                    fout.write(line)
                    continue
                fout.write(_format_atom_line(parts, parts[0], chain_map))
            elif line.startswith("TER"):
                parts = line.split()
                if len(parts) < 5:
                    fout.write(line)
                    continue
                fout.write(_format_ter_line(parts, chain_map))
            else:
                fout.write(line)
