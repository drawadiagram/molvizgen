#!/usr/bin/env python3
"""
Atom-name-based ligand selection, for figures that need to split a single
hetero-residue into sub-groups (e.g. "hot spot" atoms targeted for binder
design vs. the rest of the molecule) rather than selecting whole chains.

This complements chain-id selection (the domain/peptide convention used by
pdz_figure.py / generate_figure.py): here the unit of selection is a set of
PyMOL atom names within one ligand residue, and the source of truth for
which atoms belong to which set is an RFDiffusion3 binder-design spec JSON
(see ImpressBasePipeline's `*_binder_design.json` convention) rather than a
PDB chain id.

A binder-design spec looks like:
    {"buried":  {"input": "./input_pdbs/X.pdb", "ligand": "IAI",
                 "select_buried": {"IAI": "C1,C2,..."}},
     "partial": {"input": "./input_pdbs/X.pdb", "ligand": "IAI",
                 "select_exposed": {"IAI": "C1,..."},
                 "select_buried":  {"IAI": "C2,..."}}}
`select_buried` is always the hot-spot/targeted atom set; `select_exposed`
(only present on tasks that don't treat the whole ligand as buried) is
whatever's left of the ligand.
"""
import json

from design_spec import parse_atom_names, resolve_input_path


def atom_name_selection(base_selection, atom_names):
    """Build a PyMOL selection expression restricting `base_selection` to a
    specific set of atom names, e.g. '(complex and resn IAI) and name C1+C2+N9'.
    PyMOL's `name` selector takes a '+'-separated list."""
    if not atom_names:
        return f"none"
    return f"({base_selection}) and name {'+'.join(atom_names)}"


def load_binder_design_spec(design_json_path, task):
    """Read an RFDiffusion3 binder-design spec JSON and return
    (input_pdb_path, ligand_resn, buried_atoms, exposed_atoms) for the given
    task key ('buried' or 'partial'). `input_pdb_path` is resolved relative
    to the JSON file's directory (the convention every pipeline in this repo
    uses for its 'input' field) and made absolute. `exposed_atoms` is []
    when the task has no select_exposed block (i.e. the whole ligand is the
    buried/hot-spot target)."""
    with open(design_json_path) as f:
        spec = json.load(f)
    if task not in spec:
        raise KeyError(f"{design_json_path}: no {task!r} task (has: {list(spec)})")
    task_spec = spec[task]

    ligand = task_spec["ligand"]
    input_pdb = resolve_input_path(design_json_path, task_spec["input"])

    buried_atoms = parse_atom_names(task_spec.get("select_buried", {}).get(ligand, ""))
    exposed_atoms = parse_atom_names(task_spec.get("select_exposed", {}).get(ligand, ""))

    return input_pdb, ligand, buried_atoms, exposed_atoms
