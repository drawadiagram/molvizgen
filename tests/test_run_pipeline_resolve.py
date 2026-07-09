"""
Unit-level tests for run_pipeline.py's internal resolve() mechanics and the
figure_name/out_root/data: bookkeeping in main(), imported directly (see
conftest.py, which puts REPO_ROOT on sys.path the same way every script in
this repo reaches lib/).
"""
import os
import subprocess
import sys

import pytest
import yaml

import run_pipeline

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_resolve_anchored_data_ref_returns_root_unchanged():
    roots = {"prod": "/data/prod"}
    assert run_pipeline.resolve("${data.prod}", {}, roots) == "/data/prod"


def test_resolve_embedded_data_ref_joins_subpath():
    roots = {"prod": "/data/prod"}
    value = run_pipeline.resolve("${data.prod}/p1_in/sm_binder_design.json", {}, roots)
    assert value == "/data/prod/p1_in/sm_binder_design.json"


def test_resolve_leaves_strings_without_data_ref_untouched():
    roots = {"prod": "/data/prod"}
    assert run_pipeline.resolve("plain string", {}, roots) == "plain string"
    assert run_pipeline.resolve(42, {}, roots) == 42


def test_resolve_undeclared_data_ref_raises():
    with pytest.raises(SystemExit) as exc_info:
        run_pipeline.resolve("${data.missing}", {}, {"prod": "/data/prod"})
    assert "undeclared data root" in str(exc_info.value)


def test_resolve_step_field_refs_unaffected_by_roots_param():
    ctx = {"find": {"manifest": "/x/candidates.json", "pngs": ["/x/a.png", "/x/b.png"]}}
    assert run_pipeline.resolve("${find.manifest}", ctx, {}) == "/x/candidates.json"
    assert run_pipeline.resolve("${find.pngs}", ctx, {}) == ["/x/a.png", "/x/b.png"]
    assert run_pipeline.resolve({"in": "${find.manifest}"}, ctx, {}) == {"in": "/x/candidates.json"}


def test_resolve_unknown_step_field_ref_raises():
    with pytest.raises(SystemExit):
        run_pipeline.resolve("${missing.manifest}", {}, {})


def _run_pipeline_cli(pipeline_path, *extra_args):
    return subprocess.run(
        [sys.executable, os.path.join(REPO_ROOT, "run_pipeline.py"), str(pipeline_path), *extra_args],
        capture_output=True, text=True,
    )


def test_missing_figure_name_is_a_hard_error(tmp_path):
    pipeline_path = tmp_path / "pipeline.yaml"
    pipeline_path.write_text(yaml.dump({"steps": []}))
    result = _run_pipeline_cli(pipeline_path)
    assert result.returncode != 0
    assert "figure_name" in result.stderr


def test_step_named_data_is_rejected(tmp_path):
    pipeline = {
        "figure_name": "test_figure",
        "out_root": str(tmp_path),
        "steps": [{"name": "data", "kind": "find_flat", "args": {"dir": str(tmp_path)}}],
    }
    pipeline_path = tmp_path / "pipeline.yaml"
    pipeline_path.write_text(yaml.dump(pipeline))
    result = _run_pipeline_cli(pipeline_path)
    assert result.returncode != 0
    assert "reserved" in result.stderr


@pytest.mark.parametrize("cli_out_root,yaml_out_root,expected_subdir", [
    (None, None, "."),
    (None, "yaml_base", "yaml_base"),
    ("cli_base", "yaml_base", "cli_base"),
])
def test_out_root_precedence(tmp_path, fixtures_dir, cli_out_root, yaml_out_root, expected_subdir):
    pipeline = {
        "figure_name": "test_figure",
        "steps": [{"name": "find", "kind": "find_flat", "args": {"dir": fixtures_dir, "glob": "domain*.pdb"}}],
    }
    if yaml_out_root is not None:
        pipeline["out_root"] = str(tmp_path / yaml_out_root)
    pipeline_path = tmp_path / "pipeline.yaml"
    pipeline_path.write_text(yaml.dump(pipeline))

    extra_args = []
    if cli_out_root is not None:
        extra_args = ["--out-root", str(tmp_path / cli_out_root)]

    result = subprocess.run(
        [sys.executable, os.path.join(REPO_ROOT, "run_pipeline.py"), str(pipeline_path), *extra_args],
        capture_output=True, text=True, cwd=tmp_path,
    )
    assert result.returncode == 0, result.stderr
    expected = (tmp_path / expected_subdir / "test_figure" / "find" / "candidates.json")
    assert expected.exists()
