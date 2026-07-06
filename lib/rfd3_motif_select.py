#!/usr/bin/env python3
"""
RFD3 structured design-spec selection — a third selection convention
alongside chain ids (the domain/peptide convention used by pdz_figure.py /
generate_figure.py) and the buried/exposed ligand atom-name split
(lib/ligand_select.py, for `*_binder_design.json` specs). This one is for
the per-model design spec used by the discontinuous-scaffolds RFDiffusion3
benchmark (see that project's `scripts/mcsa_41-N.json`), shaped as:

    {"M0110_1c0p": {"input": "./mcsa_41/M0110_1c0p.pdb",
                     "ligand": "DAL,PER,FAD",
                     "contig": "A1054,280,A1335,3,A1339",
                     "select_fixed_atoms": {"A1054": "N,CA",
                                             "A1335": "O,C,CA",
                                             "A1339": "N,CA",
                                             "DAL": "N,CA,CB,C,O", ...}}}

`select_fixed_atoms` mixes two kinds of keys in one dict: protein motif
residues ("A1054" = chain "A" + PDB residue number 1054) and ligand
atom-name splits (a three-letter resn matching one of the comma-separated
tokens in "ligand"). `split_fixed_atoms` separates the two.

`contig_to_chai_positions` walks the RFD3 contig string to map each protein
motif residue to the 1-indexed sequence position it occupies in any
RFDiffusion3/Chai-1 output generated from this spec. That mapping also
holds for a *redesign* descendant of this spec (see the discontinuous-
scaffolds project's `create_redesign.py`/adaptive-branching docs): a
redesign's contig relabels which reference residues anchor the scaffold,
but preserves every free-residue run-length between them, so the
cumulative sequence-position arithmetic for a given original motif residue
comes out the same regardless of how many redesign generations produced
the fold being examined.
"""
import json
import os
import re

_PROTEIN_RES_RE = re.compile(r"^([A-Za-z])(\d+)$")


def load_design_spec(design_json_path, model_name):
    """Read an RFD3 structured design-spec JSON and return the entry for
    `model_name`, with `input` resolved to an absolute path relative to the
    JSON file's own directory (the convention every pipeline in this repo
    uses for its 'input' field)."""
    with open(design_json_path) as f:
        spec = json.load(f)
    if model_name not in spec:
        raise KeyError(f"{design_json_path}: no {model_name!r} entry (has: {list(spec)})")
    entry = dict(spec[model_name])
    base_dir = os.path.dirname(os.path.abspath(design_json_path))
    entry["input"] = os.path.normpath(os.path.join(base_dir, entry["input"]))
    return entry


def parse_atom_names(csv_str):
    """'C22,C23,N13' -> ['C22', 'C23', 'N13']. Empty/missing string -> []."""
    if not csv_str:
        return []
    return [tok.strip() for tok in csv_str.split(",") if tok.strip()]


def split_fixed_atoms(select_fixed_atoms, ligand_resns):
    """Split a select_fixed_atoms dict into (protein_atoms, ligand_atoms):
    protein_atoms: {(chain, resnum): [atom names]}, for keys matching the
    "{Chain}{ResNum}" pattern; ligand_atoms: {resn: [atom names]}, for keys
    matching one of `ligand_resns`. Raises on a key matching neither."""
    ligand_resns = set(ligand_resns)
    protein_atoms, ligand_atoms = {}, {}
    for key, atoms_csv in select_fixed_atoms.items():
        m = _PROTEIN_RES_RE.match(key)
        if key in ligand_resns:
            ligand_atoms[key] = parse_atom_names(atoms_csv)
        elif m:
            protein_atoms[(m.group(1), int(m.group(2)))] = parse_atom_names(atoms_csv)
        else:
            raise ValueError(
                f"select_fixed_atoms key {key!r} matches neither a protein "
                f"residue (Chain+ResNum) nor a declared ligand {sorted(ligand_resns)}"
            )
    return protein_atoms, ligand_atoms


def contig_to_chai_positions(contig_str):
    """Walk an RFD3 contig string and return {(chain, resnum): chai_seq_pos}
    (1-indexed): the position each named motif residue occupies in any fold
    output generated from this contig (see module docstring for why this
    stays valid across redesign generations too)."""
    mapping = {}
    pos = 1
    for token in contig_str.split(","):
        token = token.strip()
        if not token:
            continue
        m = _PROTEIN_RES_RE.match(token)
        if m:
            mapping[(m.group(1), int(m.group(2)))] = pos
            pos += 1
        else:
            pos += int(token)
    return mapping


def ligand_resn_selection(base_selection, ligand_resns):
    """Whole-ligand selection by resn, e.g. '(reference) and resn DAL+PER+FAD'."""
    resns = [r for r in ligand_resns if r]
    if not resns:
        return "none"
    return f"({base_selection}) and resn {'+'.join(resns)}"
