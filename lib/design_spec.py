#!/usr/bin/env python3
"""
Small helpers shared by the two RFDiffusion3 design-spec selection
conventions (lib/ligand_select.py's buried/exposed atom-name split, and
lib/rfd3_motif_select.py's structured per-model spec): parsing a
comma-separated atom-name list, and resolving a spec's relative "input"
path against the spec JSON's own directory (the convention every pipeline
in this repo uses for that field).
"""
import os


def parse_atom_names(csv_str):
    """'C22,C23,N13' -> ['C22', 'C23', 'N13']. Empty/missing string -> []."""
    if not csv_str:
        return []
    return [tok.strip() for tok in csv_str.split(",") if tok.strip()]


def resolve_input_path(design_json_path, relative_input):
    """Resolve a spec's relative 'input' path against the spec JSON's own
    directory, and normalize it, returning an absolute path."""
    base_dir = os.path.dirname(os.path.abspath(design_json_path))
    return os.path.normpath(os.path.join(base_dir, relative_input))
