import pytest

from pdb_normalize import detect_chain_tokens, normalize_pdb_chains


def test_detect_chain_tokens_finds_multichar_tokens_in_first_seen_order(boltz_chain_tokens_pdb):
    assert detect_chain_tokens(boltz_chain_tokens_pdb) == ["pdz", "pep"]


def test_detect_chain_tokens_on_standard_pdb_returns_single_letters(complex_ab_pdb):
    tokens = detect_chain_tokens(complex_ab_pdb)
    assert tokens == ["A", "B"]


def test_normalize_pdb_chains_rewrites_to_fixed_column_single_char(tmp_path, boltz_chain_tokens_pdb):
    out_path = tmp_path / "normalized.pdb"
    normalize_pdb_chains(boltz_chain_tokens_pdb, str(out_path), {"pdz": "A", "pep": "B"})

    lines = out_path.read_text().splitlines()
    atom_lines = [l for l in lines if l.startswith("ATOM")]
    assert len(atom_lines) == 5
    # fixed column 22 (0-indexed 21) holds the single-char chain id
    assert {l[21] for l in atom_lines[:3]} == {"A"}
    assert {l[21] for l in atom_lines[3:]} == {"B"}
    # standard fixed-column line length (up to element symbol)
    for l in atom_lines:
        assert len(l.rstrip("\n")) >= 66

    ter_lines = [l for l in lines if l.startswith("TER")]
    assert len(ter_lines) == 2
    assert ter_lines[0][21] == "A"
    assert ter_lines[1][21] == "B"


def test_normalize_pdb_chains_raises_on_unmapped_token(tmp_path, boltz_chain_tokens_pdb):
    out_path = tmp_path / "normalized.pdb"
    with pytest.raises(KeyError):
        normalize_pdb_chains(boltz_chain_tokens_pdb, str(out_path), {"pdz": "A"})


def test_normalize_pdb_chains_creates_parent_dirs(tmp_path, boltz_chain_tokens_pdb):
    out_path = tmp_path / "nested" / "dir" / "normalized.pdb"
    normalize_pdb_chains(boltz_chain_tokens_pdb, str(out_path), {"pdz": "A", "pep": "B"})
    assert out_path.exists()
