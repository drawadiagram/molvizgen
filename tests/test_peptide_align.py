import numpy as np
import pytest

from geometry import long_axis_pca
from peptide_align import backbone_coords, chain_resi_sorted, cterminal_resi, fit_cterminal_backbone, load_and_align_pair
from pymol_scene import cmd


def test_chain_resi_sorted_and_cterminal_resi(complex_ab_pdb):
    obj = "peptide_align_test_resi"
    cmd.load(complex_ab_pdb, obj)

    resi = chain_resi_sorted(obj, "B")
    assert resi == sorted(resi)
    assert len(resi) == 10  # complex_ab.pdb's chain B peptide, waters excluded

    assert cterminal_resi(obj, "B", 4) == resi[-4:]

    cmd.delete(obj)


def test_cterminal_resi_too_few_residues_exits(complex_ab_pdb):
    obj = "peptide_align_test_toofew"
    cmd.load(complex_ab_pdb, obj)

    with pytest.raises(SystemExit):
        cterminal_resi(obj, "B", 999)

    cmd.delete(obj)


def test_backbone_coords_missing_atom_exits(complex_ab_pdb):
    obj = "peptide_align_test_missing"
    cmd.load(complex_ab_pdb, obj)

    with pytest.raises(SystemExit):
        backbone_coords(obj, "Z", [1, 2, 3, 4])

    cmd.delete(obj)


def test_fit_cterminal_backbone_aligns_close_structures(complex_ab_pdb, complex_cd_pdb):
    mobile_obj = "peptide_align_test_mobile"
    target_obj = "peptide_align_test_target"
    cmd.load(complex_cd_pdb, mobile_obj)
    cmd.load(complex_ab_pdb, target_obj)

    R, t = fit_cterminal_backbone(mobile_obj, "D", target_obj, "B", 4)

    mobile_pts = backbone_coords(mobile_obj, "D", cterminal_resi(mobile_obj, "D", 4))
    target_pts = backbone_coords(target_obj, "B", cterminal_resi(target_obj, "B", 4))
    fitted = mobile_pts @ R.T + t

    # complex_ab.pdb/complex_cd.pdb are close-but-not-identical independent
    # copies of the same complex, so the fit should bring their peptides'
    # last-4-residue backbones very close (not exactly zero) together.
    rmsd = np.sqrt(np.mean(np.sum((fitted - target_pts) ** 2, axis=1)))
    assert rmsd < 1.0

    cmd.delete(mobile_obj)
    cmd.delete(target_obj)


def test_fit_cterminal_backbone_mismatched_lengths_exits(complex_ab_pdb):
    obj = "peptide_align_test_mismatch"
    cmd.load(complex_ab_pdb, obj)

    with pytest.raises(SystemExit):
        fit_cterminal_backbone(obj, "B", obj, "B", 999)

    cmd.delete(obj)


def test_load_and_align_pair_orients_design_and_aligns_reference_peptide(complex_ab_pdb, complex_cd_pdb):
    design_obj = "peptide_align_test_lap_design"
    reference_obj = "peptide_align_test_lap_reference"

    load_and_align_pair(
        complex_cd_pdb, complex_ab_pdb, "D", "A", "B", 4,
        design_obj=design_obj, reference_obj=reference_obj,
    )

    # design_obj should be oriented by the standard convention: its domain's
    # long axis vertical (Y).
    design_axis_coords = cmd.get_coords(f"{design_obj} and chain A and name CA and polymer")
    _, long_axis = long_axis_pca(design_axis_coords)
    assert abs(abs(long_axis[1]) - 1.0) < 1e-6

    # reference_obj's peptide C-terminal backbone should now sit very close
    # to design_obj's (complex_ab.pdb/complex_cd.pdb are close-but-not-
    # identical independent copies of the same complex).
    reference_pts = backbone_coords(reference_obj, "D", cterminal_resi(reference_obj, "D", 4))
    design_pts = backbone_coords(design_obj, "B", cterminal_resi(design_obj, "B", 4))
    rmsd = np.sqrt(np.mean(np.sum((reference_pts - design_pts) ** 2, axis=1)))
    assert rmsd < 1.0

    cmd.delete(design_obj)
    cmd.delete(reference_obj)
