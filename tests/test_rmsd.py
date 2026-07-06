import csv
import os

import pytest

from rmsd import (
    greedy_max_min_selection,
    load_ca_objects,
    pairwise_rmsd_matrix,
    require_single_chain_id,
    write_matrix_csv,
)


@pytest.fixture
def domain_fixtures(fixtures_dir):
    return [
        ("d1", os.path.join(fixtures_dir, "domain1.pdb")),
        ("d2", os.path.join(fixtures_dir, "domain2.pdb")),
        ("d3", os.path.join(fixtures_dir, "domain3_unrelated.pdb")),
    ]


# --- pure, no-PyMOL-calls function ---

def test_greedy_max_min_selection_picks_globally_most_dissimilar_pair_first():
    # 0 and 1 are the most dissimilar pair (distance 10); with n_select=2 that
    # pair must be selected regardless of order.
    matrix = [
        [0, 10, 1, 1],
        [10, 0, 1, 1],
        [1, 1, 0, 0.5],
        [1, 1, 0.5, 0],
    ]
    selected = greedy_max_min_selection(matrix, 2)
    assert set(selected) == {0, 1}


def test_greedy_max_min_selection_grows_by_farthest_point():
    matrix = [
        [0, 10, 1, 1],
        [10, 0, 1, 1],
        [1, 1, 0, 0.5],
        [1, 1, 0.5, 0],
    ]
    selected = greedy_max_min_selection(matrix, 3)
    assert len(selected) == 3
    assert {0, 1}.issubset(set(selected))


def test_greedy_max_min_selection_full_set():
    matrix = [[0, 1], [1, 0]]
    assert sorted(greedy_max_min_selection(matrix, 2)) == [0, 1]


def test_require_single_chain_id_returns_agreed_chain():
    candidates = [{"chain_domain": "A"}, {"chain_domain": "A"}]
    assert require_single_chain_id(candidates, "chain_domain") == "A"


def test_require_single_chain_id_exits_on_disagreement():
    candidates = [{"chain_domain": "A"}, {"chain_domain": "B"}]
    with pytest.raises(SystemExit):
        require_single_chain_id(candidates, "chain_domain")


def test_require_single_chain_id_exits_on_missing_field():
    candidates = [{"id": "x"}, {"id": "y"}]
    with pytest.raises(SystemExit):
        require_single_chain_id(candidates, "chain_domain")


def test_write_matrix_csv_roundtrip(tmp_path):
    names = ["a", "b"]
    matrix = [[0.0, 1.23456], [1.23456, 0.0]]
    out = tmp_path / "nested" / "matrix.csv"
    write_matrix_csv(str(out), names, matrix)

    with open(out, newline="") as f:
        rows = list(csv.reader(f))
    assert rows[0] == ["", "a", "b"]
    assert rows[1] == ["a", "0.0000", "1.2346"]
    assert rows[2] == ["b", "1.2346", "0.0000"]


# --- PyMOL-dependent functions ---

def test_load_ca_objects_filters_to_chain_and_ca(domain_fixtures):
    from pymol import cmd

    names = load_ca_objects(domain_fixtures[:1], "A")
    assert names == ["d1"]
    assert cmd.count_atoms("d1") > 0
    assert cmd.count_atoms("d1 and not name CA") == 0
    assert cmd.count_atoms("d1 and not chain A") == 0
    cmd.delete("d1")


def test_load_ca_objects_exits_when_chain_has_no_atoms(fixtures_dir):
    with pytest.raises(SystemExit):
        load_ca_objects([("missing_chain", os.path.join(fixtures_dir, "domain1.pdb"))], "Z")


def test_pairwise_rmsd_matrix_ordering_matches_structural_similarity(domain_fixtures):
    names = load_ca_objects(domain_fixtures, "A")
    matrix = pairwise_rmsd_matrix(names)

    assert matrix[0][0] == 0.0
    # d1 and d2 are two independently-refined copies of the same PDZ domain:
    # small but nonzero RMSD.
    assert 0 < matrix[0][1] < 1.0
    # d3 is an unrelated fold: much larger RMSD from both.
    assert matrix[0][2] > 5.0
    assert matrix[1][2] > 5.0
    # symmetric
    assert matrix[0][1] == matrix[1][0]
    assert matrix[0][2] == matrix[2][0]

    for name in names:
        from pymol import cmd
        cmd.delete(name)


def test_pairwise_rmsd_matrix_progress_callback_invoked_on_completion(domain_fixtures):
    names = load_ca_objects(domain_fixtures, "A")
    calls = []
    pairwise_rmsd_matrix(names, progress=lambda done, total: calls.append((done, total)))
    assert calls[-1] == (3, 3)  # 3 pairs among 3 structures

    for name in names:
        from pymol import cmd
        cmd.delete(name)
