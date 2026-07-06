"""
CLI-level smoke tests, one representative script per pipeline stage (plus
run_pipeline.py end to end), run exactly as a user would via subprocess.
These exercise the actual argparse surface and script-to-script contract
(manifest JSON, --flag names), not just the underlying library functions.
"""
import json
import os
import subprocess
import sys

import pytest
import yaml
from PIL import Image

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def run_script(name, *args):
    cmd = [sys.executable, os.path.join(REPO_ROOT, name), *args]
    return subprocess.run(cmd, capture_output=True, text=True)


def test_find_structures_flat_over_fixture_dir(tmp_path, fixtures_dir):
    out = tmp_path / "candidates.json"
    result = run_script("find_structures_flat.py", "--dir", fixtures_dir, "--glob", "domain*.pdb", "--out", str(out))
    assert result.returncode == 0, result.stderr

    manifest = json.loads(out.read_text())
    ids = {c["id"] for c in manifest["candidates"]}
    assert ids == {"domain1", "domain2", "domain3_unrelated"}
    assert all(os.path.isabs(c["pdb_path"]) for c in manifest["candidates"])
    assert all(c["chain_domain"] == "A" for c in manifest["candidates"])


def test_filter_top_n_over_synthetic_manifest(tmp_path):
    manifest_in = tmp_path / "in.json"
    manifest_in.write_text(json.dumps({"candidates": [
        {"id": "a", "pdb_path": "/x/a.pdb", "score": 1.0},
        {"id": "b", "pdb_path": "/x/b.pdb", "score": 3.0},
        {"id": "c", "pdb_path": "/x/c.pdb", "score": 2.0},
    ]}))
    out = tmp_path / "top.json"

    result = run_script("filter_top_n.py", "--in", str(manifest_in), "--out", str(out), "--score-field", "score", "--n", "2")
    assert result.returncode == 0, result.stderr

    kept = [c["id"] for c in json.loads(out.read_text())["candidates"]]
    assert kept == ["b", "c"]


def test_filter_diversity_over_fixture_domains(tmp_path, fixtures_dir):
    manifest_in = tmp_path / "in.json"
    manifest_in.write_text(json.dumps({"candidates": [
        {"id": "d1", "pdb_path": os.path.join(fixtures_dir, "domain1.pdb"), "chain_domain": "A"},
        {"id": "d2", "pdb_path": os.path.join(fixtures_dir, "domain2.pdb"), "chain_domain": "A"},
        {"id": "d3", "pdb_path": os.path.join(fixtures_dir, "domain3_unrelated.pdb"), "chain_domain": "A"},
    ]}))
    out_dir = tmp_path / "analysis"

    result = run_script(
        "filter_diversity.py", "--in", str(manifest_in), "--out-dir", str(out_dir),
        "--chain-field", "chain_domain", "--n-select", "2",
    )
    assert result.returncode == 0, result.stderr
    assert (out_dir / "rmsd_matrix.csv").exists()
    selection = json.loads((out_dir / "selection.json").read_text())
    assert len(selection["candidates"]) == 2


def test_pdz_figure_renders_a_nonempty_png(tmp_path, complex_ab_pdb):
    out_png = tmp_path / "out.png"
    result = run_script("pdz_figure.py", complex_ab_pdb, str(out_png), "--width", "200", "--height", "250")
    assert result.returncode == 0, result.stderr
    assert out_png.exists()
    img = Image.open(out_png)
    assert img.size == (200, 250)


def test_aligned_pair_figure_renders_two_nonempty_pngs(tmp_path, complex_ab_pdb, complex_cd_pdb):
    out_reference = tmp_path / "reference.png"
    out_design = tmp_path / "design.png"
    result = run_script(
        "aligned_pair_figure.py", complex_cd_pdb, complex_ab_pdb, str(out_reference), str(out_design),
        "--reference-chain-domain", "C", "--reference-chain-peptide", "D",
        "--design-chain-domain", "A", "--design-chain-peptide", "B",
        "--width", "200", "--height", "250",
    )
    assert result.returncode == 0, result.stderr
    assert out_reference.exists()
    assert out_design.exists()
    assert Image.open(out_reference).size == (200, 250)
    assert Image.open(out_design).size == (200, 250)


def test_aligned_overlay_figure_renders_a_nonempty_png(tmp_path, complex_ab_pdb, complex_cd_pdb):
    out_png = tmp_path / "overlay.png"
    result = run_script(
        "aligned_overlay_figure.py", complex_cd_pdb, complex_ab_pdb, str(out_png),
        "--reference-chain-domain", "C", "--reference-chain-peptide", "D",
        "--design-chain-domain", "A", "--design-chain-peptide", "B",
        "--width", "200", "--height", "250",
    )
    assert result.returncode == 0, result.stderr
    assert out_png.exists()
    assert Image.open(out_png).size == (200, 250)


def test_aligned_overlay_figure_also_writes_solo_views(tmp_path, complex_ab_pdb, complex_cd_pdb):
    out_png = tmp_path / "overlay.png"
    out_reference = tmp_path / "reference_alone.png"
    out_design = tmp_path / "design_alone.png"
    result = run_script(
        "aligned_overlay_figure.py", complex_cd_pdb, complex_ab_pdb, str(out_png),
        "--reference-chain-domain", "C", "--reference-chain-peptide", "D",
        "--design-chain-domain", "A", "--design-chain-peptide", "B",
        "--out-reference", str(out_reference), "--out-design", str(out_design),
        "--width", "200", "--height", "250",
    )
    assert result.returncode == 0, result.stderr
    assert out_png.exists()
    assert out_reference.exists()
    assert out_design.exists()


def test_montage_figures_tiles_two_pngs(tmp_path):
    a = tmp_path / "a.png"
    b = tmp_path / "b.png"
    Image.new("RGB", (50, 50), "red").save(a)
    Image.new("RGB", (50, 50), "blue").save(b)
    out = tmp_path / "montage.png"

    result = run_script("montage_figures.py", str(a), str(b), "--rows", "1", "--cols", "2", "--out", str(out), "--no-scale", "--no-trim")
    assert result.returncode == 0, result.stderr
    assert out.exists()


def test_run_pipeline_end_to_end(tmp_path, fixtures_dir):
    pipeline = {
        "out_dir": str(tmp_path / "out"),
        "steps": [
            {"name": "find", "kind": "find_flat", "args": {"dir": fixtures_dir, "glob": "domain*.pdb"}},
            {"name": "select", "kind": "filter_diversity", "args": {
                "in": "${find.manifest}", "chain_field": "chain_domain", "n_select": 2,
            }},
            {"name": "heatmap", "kind": "plot_heatmap", "args": {
                "matrix": "${select.matrix}", "out": "heatmap.png",
            }},
        ],
    }
    pipeline_path = tmp_path / "pipeline.yaml"
    pipeline_path.write_text(yaml.dump(pipeline))

    result = subprocess.run(
        [sys.executable, os.path.join(REPO_ROOT, "run_pipeline.py"), str(pipeline_path)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "out" / "heatmap" / "heatmap.png").exists()
    assert (tmp_path / "out" / "select" / "rmsd_matrix.csv").exists()
