#!/usr/bin/env python3
"""
Fourth selection/alignment convention, alongside chain-id whole-chain
selection, lib/ligand_select.py's atom-name ligand split, and
lib/rfd3_motif_select.py's contig-based motif-residue mapping: rigid-body
Kabsch superposition (lib/kabsch.py) of one peptide chain's C-terminal
residues onto another peptide chain's C-terminal residues, matched purely by
position along the chain (N -> C order) rather than residue identity or an
external contig map.

This is the right tool when a short peptide fragment (e.g. a PDZ domain's
minimal conserved C-terminal binding motif, used as an RFDiffusion/design
template) and a longer peptide produced downstream by the same pipeline are
known to share that same binding motif at their C-terminus, so the last N
residues of each correspond 1:1 in N -> C order — no sequence/structure
alignment step (lib/rmsd.py's cealign) or external correspondence map
(lib/rfd3_motif_select.py's contig) is needed to find the pairing.

`load_and_align_pair` bundles the full recipe both figure scripts built on
this module need before they diverge on how to render it
(aligned_pair_figure.py: two separate panels; aligned_overlay_figure.py: one
combined panel) — load `design`, orient it with the standard
long-axis-vertical convention, load `reference`, then Kabsch-fit it into
`design`'s now-oriented frame.
"""
import sys

import numpy as np

from kabsch import as_homogeneous, kabsch_fit
from orient import orient_long_axis_vertical
from pymol_scene import cmd

BACKBONE_ATOMS = ("N", "CA", "C", "O")


def chain_resi_sorted(obj, chain):
    """Every polymer residue number present for `chain` in `obj`, ascending
    (N -> C-terminus under standard PDB numbering)."""
    model = cmd.get_model(f"{obj} and chain {chain} and polymer and name CA")
    return sorted({int(a.resi) for a in model.atom})


def cterminal_resi(obj, chain, n_residues):
    """The last `n_residues` residue numbers of `chain` in `obj`, in N -> C
    order — e.g. a PDZ-binding motif's conserved C-terminal residues,
    regardless of how many residues precede them or what numbering offset
    the structure uses."""
    resi = chain_resi_sorted(obj, chain)
    if len(resi) < n_residues:
        sys.exit(f"{obj!r} chain {chain} has only {len(resi)} residue(s), need >= {n_residues}")
    return resi[-n_residues:]


def backbone_coords(obj, chain, resi_list, atom_names=BACKBONE_ATOMS):
    """Ordered Nx3 backbone coordinates for `resi_list` (ints) of `chain` in
    `obj`, one row per (resi, atom_name) pair in iteration order. Exits
    (rather than falling back, as lib/rfd3_motif_select.py's hot-spot atoms
    do) if a named atom isn't present exactly once, since a peptide backbone
    is expected to be fully resolved."""
    pts = []
    for resi in resi_list:
        for atom_name in atom_names:
            sel = f"{obj} and chain {chain} and resi {resi} and name {atom_name}"
            n = cmd.count_atoms(sel)
            if n != 1:
                sys.exit(f"{sel!r} matched {n} atom(s) (need exactly 1) — can't build a backbone Kabsch fit")
            pts.append(cmd.get_coords(sel)[0])
    return np.array(pts)


def fit_cterminal_backbone(mobile_obj, mobile_chain, target_obj, target_chain, n_residues):
    """(R, t) Kabsch-superposing `mobile_obj`'s last `n_residues` peptide
    backbone atoms onto `target_obj`'s, matched by position (N -> C order)
    rather than residue identity or an external contig map — apply with
    `cmd.transform_selection` (via lib/kabsch.py's `as_homogeneous`) to
    `mobile_obj` as a whole to bring its entire structure into
    `target_obj`'s coordinate frame."""
    mobile_resi = cterminal_resi(mobile_obj, mobile_chain, n_residues)
    target_resi = cterminal_resi(target_obj, target_chain, n_residues)
    mobile_pts = backbone_coords(mobile_obj, mobile_chain, mobile_resi)
    target_pts = backbone_coords(target_obj, target_chain, target_resi)
    return kabsch_fit(mobile_pts, target_pts)


def load_and_align_pair(reference_pdb, design_pdb, reference_peptide_chain,
                         design_domain_chain, design_peptide_chain,
                         align_length, design_obj="design", reference_obj="reference"):
    """Load `design_pdb`/`reference_pdb` into the current PyMOL session as
    `design_obj`/`reference_obj`, orient `design_obj` with the standard
    domain/peptide long-axis-vertical convention (lib/orient.py), then
    rigidly Kabsch-fit `reference_obj` (its whole structure, both chains)
    into `design_obj`'s now-oriented frame via `fit_cterminal_backbone`. No
    `reference_domain_chain` is needed: the alignment only touches
    `reference_peptide_chain`, and it carries the rest of `reference_obj`
    (including its domain chain, whatever id it uses) along rigidly."""
    cmd.load(design_pdb, design_obj)
    cmd.remove(f"{design_obj} and hydro")
    cmd.dss(design_obj)

    orient_long_axis_vertical(
        f"{design_obj} and chain {design_domain_chain} and name CA and polymer",
        f"{design_obj} and chain {design_peptide_chain} and name CA and polymer",
        design_obj,
    )

    cmd.load(reference_pdb, reference_obj)
    cmd.remove(f"{reference_obj} and hydro")
    cmd.dss(reference_obj)

    R, t = fit_cterminal_backbone(
        reference_obj, reference_peptide_chain, design_obj, design_peptide_chain, align_length,
    )
    cmd.transform_selection(reference_obj, as_homogeneous(R, t), homogenous=1)
