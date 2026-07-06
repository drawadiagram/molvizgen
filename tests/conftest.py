import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LIB_DIR = os.path.join(REPO_ROOT, "lib")
FIXTURES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")

# Every script in this repo reaches lib/ modules via this same path-hack
# (see CLAUDE.md); tests import lib/ modules the same way scripts do.
if LIB_DIR not in sys.path:
    sys.path.insert(0, LIB_DIR)
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import pytest


@pytest.fixture(scope="session")
def repo_root():
    return REPO_ROOT


@pytest.fixture(scope="session")
def fixtures_dir():
    return FIXTURES_DIR


@pytest.fixture(scope="session")
def complex_ab_pdb():
    """Real PDZ-domain (chain A) + peptide (chain B) complex."""
    return os.path.join(FIXTURES_DIR, "complex_ab.pdb")


@pytest.fixture(scope="session")
def complex_cd_pdb():
    """A second, independently-refined copy of the same PDZ-peptide complex
    (chain C domain + chain D peptide) - structurally very close to
    complex_ab.pdb but not identical, for small-nonzero-RMSD comparisons."""
    return os.path.join(FIXTURES_DIR, "complex_cd.pdb")


@pytest.fixture(scope="session")
def unrelated_lipase_pdb():
    """An unrelated protein fold (a lipase fragment, chain A), for
    large-RMSD/dissimilarity comparisons against the PDZ domain fixtures."""
    return os.path.join(FIXTURES_DIR, "unrelated_lipase.pdb")


@pytest.fixture(scope="session")
def domain1_pdb():
    """Single-chain-A-only PDZ domain, CA+backbone, from complex_ab.pdb."""
    return os.path.join(FIXTURES_DIR, "domain1.pdb")


@pytest.fixture(scope="session")
def domain2_pdb():
    """A second, independently-refined copy of the same PDZ domain
    (relabeled to chain A), close-but-not-identical to domain1.pdb."""
    return os.path.join(FIXTURES_DIR, "domain2.pdb")


@pytest.fixture(scope="session")
def domain3_unrelated_pdb():
    """An unrelated fold (lipase fragment, chain A), far in RMSD from
    domain1.pdb/domain2.pdb."""
    return os.path.join(FIXTURES_DIR, "domain3_unrelated.pdb")


@pytest.fixture(scope="session")
def boltz_chain_tokens_pdb():
    """A hand-crafted PDB with Boltz-style multi-character whitespace-
    tokenized chain ids ('pdz'/'pep') instead of fixed-column single-char
    chain ids."""
    return os.path.join(FIXTURES_DIR, "boltz_chain_tokens.pdb")


@pytest.fixture(scope="session")
def binder_design_spec_json():
    return os.path.join(FIXTURES_DIR, "binder_design_spec.json")


@pytest.fixture(scope="session")
def rfd3_design_spec_json():
    return os.path.join(FIXTURES_DIR, "rfd3_design_spec.json")
