import json
import sys

import pytest

from manifest import read_manifest, read_manifest_stdin_or_path, write_manifest


def test_write_then_read_roundtrip(tmp_path):
    candidates = [
        {"id": "a", "pdb_path": "/abs/a.pdb", "score": 1.5},
        {"id": "b", "pdb_path": "/abs/b.pdb", "score": 2.5},
    ]
    out = tmp_path / "nested" / "manifest.json"
    write_manifest(str(out), candidates)

    assert out.exists()
    on_disk = json.loads(out.read_text())
    assert on_disk == {"candidates": candidates}
    assert read_manifest(str(out)) == candidates


def test_read_manifest_stdin(monkeypatch):
    candidates = [{"id": "x", "pdb_path": "/x.pdb"}]
    monkeypatch.setattr(sys, "stdin", __import__("io").StringIO(json.dumps({"candidates": candidates})))
    assert read_manifest_stdin_or_path("-") == candidates


def test_read_manifest_stdin_or_path_delegates_to_file(tmp_path):
    candidates = [{"id": "y", "pdb_path": "/y.pdb"}]
    path = tmp_path / "m.json"
    write_manifest(str(path), candidates)
    assert read_manifest_stdin_or_path(str(path)) == candidates


def test_read_manifest_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_manifest(str(tmp_path / "does_not_exist.json"))
