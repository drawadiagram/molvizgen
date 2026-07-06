#!/usr/bin/env python3
"""
Point-correspondence rigid-body superposition (Kabsch algorithm).

This complements lib/rmsd.py's whole-chain, cealign-based structural-
diversity core: cealign finds its own residue correspondence from a
sequence/structure alignment, which is the right tool when two structures
are "the same kind of thing" solved/predicted independently. Some figures
instead already know the correspondence explicitly — e.g. a folded design's
sequence position N is understood (via an RFDiffusion3 contig) to be the
same residue as reference-structure residue (chain, resnum) — and just need
the least-squares rotation/translation that best superimposes one paired
point set onto another. That's what this module provides.
"""
import numpy as np


def kabsch_fit(mobile, target):
    """Return (R, t) such that `R @ mobile[i] + t` best approximates
    `target[i]` in the least-squares sense, for paired Nx3 point arrays
    `mobile`/`target` (N >= 3, same order/length)."""
    mobile = np.asarray(mobile, dtype=float)
    target = np.asarray(target, dtype=float)
    if mobile.shape != target.shape or len(mobile) < 3:
        raise ValueError(
            f"kabsch_fit needs >=3 paired points of matching shape, "
            f"got mobile={mobile.shape} target={target.shape}"
        )

    mobile_c = mobile.mean(axis=0)
    target_c = target.mean(axis=0)
    H = (mobile - mobile_c).T @ (target - target_c)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T)) or 1.0
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    t = target_c - R @ mobile_c
    return R, t


def as_homogeneous(R, t):
    """Pack (R, t) into the flattened row-major 4x4 matrix
    `cmd.transform_selection(..., homogenous=1)` expects."""
    m4 = np.eye(4)
    m4[:3, :3] = R
    m4[:3, 3] = t
    return list(m4.flatten())
