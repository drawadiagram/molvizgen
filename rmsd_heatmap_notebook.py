import marimo

__generated_with = "0.19.11"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _():
    import json
    import subprocess
    import sys
    from pathlib import Path

    return Path, json, subprocess, sys


@app.cell
def _(Path):
    REPO_DIR = Path(__file__).resolve().parent
    DEFAULT_DATA_DIR = REPO_DIR.parent / "new_decoys"
    return DEFAULT_DATA_DIR, REPO_DIR


@app.cell
def _(mo):
    mo.md(r"""
    # Pairwise RMSD diversity pipeline

    This notebook documents the YAML pipeline format consumed by
    `run_pipeline.py`, then runs a real instance of that pipeline end to
    end against an external structure directory: FIND every structure,
    FILTER it down to a pairwise-dissimilar subset by chain-CA RMSD
    (`lib/rmsd.py`), and PLOT the full pairwise RMSD matrix as a heatmap.

    molvizgen is a tool repo, not a data repo — it takes a directory of PDB
    files as input (below, defaulting to the sibling `../new_decoys`
    directory). The three CLI scripts involved — `find_structures_flat.py`,
    `filter_diversity.py`, `figures/plot_rmsd_heatmap.py` — are unchanged by
    this notebook; it drives them exactly as `run_pipeline.py` would, and
    displays their output inline.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## The YAML pipeline format

    A pipeline file is a top-level mapping with two keys:

    | Key | Meaning |
    |---|---|
    | `out_dir` | Base directory for step outputs (default `.`). Each step gets its own subdirectory, `out_dir/<step name>/`. |
    | `steps` | A list of steps, run in order. |

    Each step is:

    ```yaml
    - name: <unique step name>
      kind: <one of the step kinds below>
      args:
        <key>: <value>          # mapped 1:1 to --<key-with-dashes> on the underlying script
    ```

    A later step can reference an earlier step's declared output with
    `${step_name.field}` — e.g. `${find.manifest}`. The reference must be
    the *entire* string value (not embedded in a larger string); it
    resolves to whatever type that field holds, a string path or a list
    of paths.

    ### Step kinds

    | `kind` | Key args | Declares |
    |---|---|---|
    | `find_flat` | `dir`, `glob` (default `*.pdb`), `chain_domain`, `chain_peptide` | `manifest` |
    | `find_campaign` | `prod_root`, `groups`, `raw_chain_domain`, `raw_chain_peptide`, `normalized_dir` | `manifest` |
    | `filter_best_score` | `in`, `key`, `score_field`, `mode` (`max`\|`min`) | `manifest` |
    | `filter_diversity` | `in`, `chain_field` (default `chain_domain`), `n_select` | `manifest`, `out_dir`, `matrix` (the `rmsd_matrix.csv` this step wrote) |
    | `plot_heatmap` | `matrix` (a `rmsd_matrix.csv`) **or** `in` + `chain_field` (compute from a manifest), `title`, `dpi`, `annotate` | `image` |
    | `generate_each` | `selection`, `script` (default `figures/pdz_figure.py`), plus any extra flags forwarded per-candidate | `pngs`, `out_dir` |
    | `assemble` | `images` (a list, or a single `${...}` ref to one), `rows`, `cols`, `out`, `target_width_in`, `dpi`, `padding`, `bg`, `no_scale` | `image` |

    `filter_diversity` is a thin wrapper around `lib/rmsd.py`: it loads
    chain-`chain_field` CA atoms for every candidate
    (`load_ca_objects`), computes the full symmetric matrix
    (`pairwise_rmsd_matrix`, via PyMOL `cealign`), writes it to
    `rmsd_matrix.csv`, then runs the greedy max-min selection
    (`greedy_max_min_selection`) to pick an `n_select`-sized
    pairwise-dissimilar subset. `plot_heatmap` renders that same matrix
    (or computes an equivalent one directly from a manifest) as a PNG.
    """)
    return


@app.cell
def _(REPO_DIR, mo):
    _yaml_text = (REPO_DIR / "rmsd_heatmap_pipeline.yaml").read_text()
    mo.md(
        f"""
        ### A concrete example

        `rmsd_heatmap_pipeline.yaml` — FIND every `*.pdb` in the sibling
        `../new_decoys` directory, FILTER to 5 pairwise-dissimilar
        structures, PLOT the full matrix as a heatmap. This is the same
        find → filter → plot sequence the rest of this notebook runs
        directly (not via `run_pipeline.py`, so it can show progress and
        render results inline).

        ```yaml
        {_yaml_text}
        ```
        """
    )
    return


@app.cell
def _(mo):
    mo.md("""
    ## Run it
    """)
    return


@app.cell
def _(DEFAULT_DATA_DIR, mo):
    pdb_dir_input = mo.ui.text(value=str(DEFAULT_DATA_DIR), label="PDB directory", full_width=True)
    n_select_input = mo.ui.number(value=5, start=2, stop=20, step=1, label="n_select")
    run_button = mo.ui.run_button(label="Run find → filter_diversity → plot_heatmap")
    mo.vstack([pdb_dir_input, mo.hstack([n_select_input, run_button], justify="start", gap=2)])
    return n_select_input, pdb_dir_input, run_button


@app.cell
def _(
    REPO_DIR,
    mo,
    n_select_input,
    pdb_dir_input,
    run_button,
    subprocess,
    sys,
):
    mo.stop(not run_button.value, mo.md("*Click **Run** above to execute the pipeline against the PDB directory given.*"))

    out_dir = REPO_DIR / "analysis" / "notebook_demo"
    out_dir.mkdir(parents=True, exist_ok=True)

    candidates_path = out_dir / "candidates.json"
    with mo.status.spinner(title="FIND: enumerating *.pdb files ..."):
        subprocess.run(
            [
                sys.executable, str(REPO_DIR / "find_structures_flat.py"),
                "--dir", pdb_dir_input.value,
                "--out", str(candidates_path),
            ],
            check=True, capture_output=True, text=True,
        )

    with mo.status.spinner(title="FILTER: pairwise chain-CA RMSD + greedy max-min selection ..."):
        subprocess.run(
            [
                sys.executable, str(REPO_DIR / "filter_diversity.py"),
                "--in", str(candidates_path),
                "--chain-field", "chain_domain",
                "--n-select", str(int(n_select_input.value)),
                "--out-dir", str(out_dir),
            ],
            check=True, capture_output=True, text=True,
        )

    heatmap_path = out_dir / "heatmap.png"
    with mo.status.spinner(title="PLOT: rendering the RMSD heatmap ..."):
        subprocess.run(
            [
                sys.executable, str(REPO_DIR / "figures" / "plot_rmsd_heatmap.py"),
                "--matrix", str(out_dir / "rmsd_matrix.csv"),
                "--out", str(heatmap_path),
                "--title", "Pairwise Cα RMSD",
            ],
            check=True, capture_output=True, text=True,
        )
    return heatmap_path, out_dir


@app.cell
def _(heatmap_path, json, mo, out_dir):
    _selection = json.loads((out_dir / "selection.json").read_text())
    _ids = [c["id"] for c in _selection["candidates"]]
    mo.vstack(
        [
            mo.md("### Selected pairwise-dissimilar subset\n" + "\n".join(f"- `{cid}`" for cid in _ids)),
            mo.image(str(heatmap_path)),
        ]
    )
    return


if __name__ == "__main__":
    app.run()
