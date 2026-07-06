import os

import pytest

from ligand_select import atom_name_selection, load_binder_design_spec, parse_atom_names


def test_parse_atom_names_splits_and_strips():
    assert parse_atom_names("C22,C23, N13 ") == ["C22", "C23", "N13"]


def test_parse_atom_names_empty_or_missing():
    assert parse_atom_names("") == []
    assert parse_atom_names(None) == []


def test_atom_name_selection_builds_plus_joined_name_selector():
    sel = atom_name_selection("complex and resn IAI", ["C1", "C2", "N9"])
    assert sel == "(complex and resn IAI) and name C1+C2+N9"


def test_atom_name_selection_empty_atom_names_is_none():
    assert atom_name_selection("complex", []) == "none"


def test_load_binder_design_spec_buried_task(binder_design_spec_json, fixtures_dir):
    input_pdb, ligand, buried, exposed = load_binder_design_spec(binder_design_spec_json, "buried")
    assert ligand == "IAI"
    assert buried == ["C1", "C2", "N9"]
    assert exposed == []
    assert os.path.isabs(input_pdb)
    assert input_pdb == os.path.normpath(os.path.join(fixtures_dir, "complex_ab.pdb"))


def test_load_binder_design_spec_partial_task(binder_design_spec_json):
    input_pdb, ligand, buried, exposed = load_binder_design_spec(binder_design_spec_json, "partial")
    assert buried == ["C1", "C2"]
    assert exposed == ["C3", "C4"]


def test_load_binder_design_spec_unknown_task_raises(binder_design_spec_json):
    with pytest.raises(KeyError):
        load_binder_design_spec(binder_design_spec_json, "nonexistent")
