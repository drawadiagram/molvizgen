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
import re

from design_spec import parse_atom_names, resolve_input_path

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
    entry["input"] = resolve_input_path(design_json_path, entry["input"])
    return entry


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


def parse_anchor_residues(anchor_residues_csv):
    """Parse the discontinuous-scaffolds project's analysis.py Step 8
    `anchor_residues` CSV cell (e.g. "A54,A56,A58,A59": motif residues whose
    backbone atoms all fell under the RMSD threshold, i.e. candidates a
    subsequent redesign would hold fixed) into a set of (chain, resnum)
    tuples — the same key shape `contig_to_chai_positions` and
    `split_fixed_atoms` use. Empty string -> empty set."""
    residues = set()
    for token in anchor_residues_csv.split(","):
        token = token.strip()
        if not token:
            continue
        m = _PROTEIN_RES_RE.match(token)
        if not m:
            raise ValueError(f"anchor residue token {token!r} doesn't match Chain+ResNum")
        residues.add((m.group(1), int(m.group(2))))
    return residues


def anchor_chai_positions(anchor_residues_csv, contig_str):
    """Map an `anchor_residues` cell to the set of chai sequence positions
    those residues occupy under `contig_str` (that generation's own contig).
    Because `contig_to_chai_positions`'s position arithmetic is invariant
    across redesign generations (see module docstring), this position set
    stays valid for identifying the *same* physical residues in any later
    redesign generation's fold too — even though that generation's own
    `select_fixed_atoms` keys are renumbered relative to the original
    reference PDB (an anchor residue carried over from a prior generation is
    renamed to its chai sequence position, e.g. an original "A59" surviving
    as a redesign's "A36"; a fresh, non-anchor reference residue is instead
    renumbered starting at 900 — see the discontinuous-scaffolds project's
    `create_redesign.py`). This is what lets a figure highlight "the portion
    that was subject to anchoring" in a redesign's fold using only the
    anchor residues identified in an *earlier* generation."""
    contig_map = contig_to_chai_positions(contig_str)
    return {contig_map[res] for res in parse_anchor_residues(anchor_residues_csv) if res in contig_map}
