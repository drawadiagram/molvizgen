import os

from design_spec import parse_atom_names, resolve_input_path


def test_parse_atom_names_splits_and_strips():
    assert parse_atom_names("C22,C23, N13 ") == ["C22", "C23", "N13"]


def test_parse_atom_names_empty_or_missing():
    assert parse_atom_names("") == []
    assert parse_atom_names(None) == []


def test_resolve_input_path_relative_to_json_dir(fixtures_dir):
    json_path = os.path.join(fixtures_dir, "binder_design_spec.json")
    resolved = resolve_input_path(json_path, "./complex_ab.pdb")
    assert os.path.isabs(resolved)
    assert resolved == os.path.normpath(os.path.join(fixtures_dir, "complex_ab.pdb"))


def test_resolve_input_path_handles_parent_dir_references(fixtures_dir):
    json_path = os.path.join(fixtures_dir, "sub", "spec.json")
    resolved = resolve_input_path(json_path, "../complex_ab.pdb")
    assert resolved == os.path.normpath(os.path.join(fixtures_dir, "complex_ab.pdb"))
