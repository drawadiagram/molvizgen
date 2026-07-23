#!/usr/bin/env python3
"""
cealign against a symmetric-repeat assembly, a variant of lib/rmsd.py's
cealign usage for when the two structures being compared are not the same
size: a single de novo backbone (RFDiffusion3 output, one protomer) vs. a
downstream fold prediction of the *whole* oligomeric assembly it seeds (a
production convention observed in the small_molecule_binding campaigns here:
ColabFold/AF2 predicts a homo-oligomer as one long single-chain object, N
copies of the protomer concatenated back-to-back with continuous residue
numbering, rather than N separate chains).

Plain cealign(target, mobile) over the full, mismatched-length selections
lets CE pick whatever local window of `target` happens to score best, which
is not necessarily the repeat unit that actually corresponds to `mobile` —
so this scans fixed-size windows along `target` explicitly (each tried via a
disposable `cmd.create` copy, same "probe, don't disturb the real mobile
object" idiom lib/rmsd.py's pairwise_rmsd_matrix uses) and keeps whichever
window cealign scores best (longest alignment, then lowest RMSD as
tiebreak), before applying that one alignment for real.
"""
from pymol_scene import cmd


def find_best_repeat_window(mobile_sel, target_obj, target_chain, window_size, stride=None):
    """Slide a `window_size`-residue window (by CA count, not raw resi
    arithmetic, so it's robust to gaps) along `target_obj`'s `target_chain`,
    cealign a disposable copy of `mobile_sel` onto each window, and return
    the best-scoring window as `(lo_resi, hi_resi, cealign_result_dict)`."""
    stride = stride or window_size
    model = cmd.get_model(f"{target_obj} and chain {target_chain} and name CA")
    resis = sorted({int(a.resi) for a in model.atom})

    best = None
    for start in range(0, len(resis) - window_size + 1, stride):
        lo, hi = resis[start], resis[start + window_size - 1]
        window_sel = f"{target_obj} and chain {target_chain} and resi {lo}-{hi} and name CA"

        cmd.create("__oligomer_align_probe", mobile_sel)
        try:
            result = cmd.cealign(window_sel, "__oligomer_align_probe")
        finally:
            cmd.delete("__oligomer_align_probe")

        score = (result["alignment_length"], -result["RMSD"])
        if best is None or score > best[0]:
            best = (score, lo, hi, result)

    if best is None:
        raise ValueError(f"{target_obj}/chain {target_chain}: fewer than {window_size} CA atoms, no window to try")
    _, lo, hi, result = best
    return lo, hi, result


def align_onto_best_repeat(mobile_sel, target_obj, target_chain, window_size, stride=None):
    """Find the best-matching window (`find_best_repeat_window`) and then
    actually cealign `mobile_sel`'s real object onto it, moving it into
    `target_obj`'s coordinate frame. Returns `(lo_resi, hi_resi,
    cealign_result_dict)` for the real, applied alignment."""
    lo, hi, _ = find_best_repeat_window(mobile_sel, target_obj, target_chain, window_size, stride)
    window_sel = f"{target_obj} and chain {target_chain} and resi {lo}-{hi} and name CA"
    result = cmd.cealign(window_sel, mobile_sel)
    return lo, hi, result
