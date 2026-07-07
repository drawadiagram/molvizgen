import os

import pytest

from rfd3_motif_select import (
    anchor_chai_positions,
    contig_to_chai_positions,
    ligand_resn_selection,
    load_design_spec,
    parse_anchor_residues,
    parse_atom_names,
    split_fixed_atoms,
)


def test_parse_atom_names_matches_ligand_select_behavior():
    assert parse_atom_names("N,CA,CB") == ["N", "CA", "CB"]
    assert parse_atom_names("") == []


def test_load_design_spec_resolves_input_relative_to_json_dir(rfd3_design_spec_json, fixtures_dir):
    entry = load_design_spec(rfd3_design_spec_json, "M0110_1c0p")
    assert entry["ligand"] == "DAL,PER,FAD"
    assert os.path.isabs(entry["input"])
    assert entry["input"] == os.path.normpath(os.path.join(fixtures_dir, "complex_ab.pdb"))


def test_load_design_spec_unknown_model_raises(rfd3_design_spec_json):
    with pytest.raises(KeyError):
        load_design_spec(rfd3_design_spec_json, "nonexistent_model")


def test_split_fixed_atoms_separates_protein_and_ligand_keys():
    select_fixed_atoms = {
        "A11": "N,CA",
        "A35": "O,C,CA",
        "DAL": "N,CA,CB,C,O",
    }
    protein_atoms, ligand_atoms = split_fixed_atoms(select_fixed_atoms, ["DAL", "PER", "FAD"])
    assert protein_atoms == {
        ("A", 11): ["N", "CA"],
        ("A", 35): ["O", "C", "CA"],
    }
    assert ligand_atoms == {"DAL": ["N", "CA", "CB", "C", "O"]}


def test_split_fixed_atoms_raises_on_unrecognized_key():
    with pytest.raises(ValueError):
        split_fixed_atoms({"???": "N,CA"}, ["DAL"])


def test_contig_to_chai_positions_walks_free_and_motif_runs():
    # A11, then 20 free residues, then A35, then 3 free, then A40
    mapping = contig_to_chai_positions("A11,20,A35,3,A40")
    assert mapping[("A", 11)] == 1
    assert mapping[("A", 35)] == 1 + 20 + 1
    assert mapping[("A", 40)] == 1 + 20 + 1 + 3 + 1


def test_ligand_resn_selection_builds_resn_selector():
    sel = ligand_resn_selection("reference", ["DAL", "PER", "FAD"])
    assert sel == "(reference) and resn DAL+PER+FAD"


def test_ligand_resn_selection_empty_is_none():
    assert ligand_resn_selection("reference", []) == "none"


def test_parse_anchor_residues_parses_csv_cell():
    assert parse_anchor_residues("A54,A56,A58,A59") == {("A", 54), ("A", 56), ("A", 58), ("A", 59)}


def test_parse_anchor_residues_empty_is_empty_set():
    assert parse_anchor_residues("") == set()


def test_parse_anchor_residues_raises_on_malformed_token():
    with pytest.raises(ValueError):
        parse_anchor_residues("A54,???")


def test_anchor_chai_positions_maps_through_contig():
    # Same contig as the module docstring example: A54, A56, A58, A59 sit at
    # chai sequence positions 1, 2, 3, 4 (no free runs between them).
    contig = "A54,0,A56,0,A58,0,A59"
    assert anchor_chai_positions("A59", contig) == {4}
    assert anchor_chai_positions("A54,A59", contig) == {1, 4}


def test_anchor_chai_positions_ignores_residues_absent_from_contig():
    assert anchor_chai_positions("A999", "A54,0,A56") == set()


def test_anchor_chai_positions_empty_anchor_residues_is_empty_set():
    assert anchor_chai_positions("", "A54,0,A56") == set()
