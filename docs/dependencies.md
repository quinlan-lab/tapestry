# Python dependencies

The shared Docker image supplies the Python runtime for both Nextflow workflows.
Container users do not need a local Python environment. For local work, use
Python 3.11 and a `.venv`; dependency files have these roles:

| File | Scope |
|---|---|
| [`requirements.txt`](../requirements.txt) | Shared runtime for trio and pedigree processing |
| [`requirements-ci.txt`](../requirements-ci.txt) | Pyright, installed alongside the runtime in CI |
| [`requirements-dev.txt`](../requirements-dev.txt) | Includes both files above, plus local notebook, figure, IGV-server, and legacy-script dependencies |

## Local development

Start with a fresh `.venv`; keep the previous environment separately if you
need it for older analyses. Installing a shorter requirements list into an
existing environment does not remove its other packages. From the repository root:

```bash
python3.11 -m venv .venv
.venv/bin/python -m pip install -r requirements-dev.txt
PYTHONPATH=src:src/util .venv/bin/python your_script.py
.venv/bin/pyright
```

CairoSVG also needs a system Cairo library. Install Cairo and any required
build dependencies for your platform using the
[CairoSVG installation instructions](https://cairosvg.org/documentation/#installation).

The Dockerfile builds pysam and cyvcf2 from source to use system OpenSSL on
FIPS-enabled hosts. A normal local pip installation can select binary wheels;
use the supported container for that cluster runtime behavior.

## Dependency audit

The consolidation checked imports in 82 Python files and all code cells in four
notebooks, plus shell invocations and optional file-format engines. No notebook
cells were skipped because of parsing errors. This is an audit of checked-in
repository uses, not of other projects installed in the previous environment.

- Merged `requirements-pipeline.txt` and `requirements-trio.txt` into the shared
  `requirements.txt`, preserving their 19 package versions.
- Added an explicit Matplotlib 3.11.2 pin because runtime plotting code imports
  it directly; this matches the already-tested container version.
- Kept Pyright separate and included it in the developer environment.
- Replaced the former 195-package environment capture with direct developer
  requirements and the notebook tools listed below. Transitive dependencies are
  resolved by pip; these lists are not a complete transitive lockfile. The old
  environment capture remains in Git history.

| Developer dependency | Repository use |
|---|---|
| aiohttp | `serve_data_for_igv.py` |
| CairoSVG | SVG-to-PDF conversion in `manuscript/fig1.py`; missing from the old requirements file |
| memory-profiler | Legacy all-CpG scripts and shell profiling commands |
| openpyxl | Excel inputs read through pandas in `src/util/get_id_to_paths.py` |
| Pillow | Raster-image work in `wiki/motivation/trio_discovery.py`; version 13.0.0 matches the existing container |
| Jupyter, JupyterLab, Notebook, ipykernel, nbconvert, ipywidgets | Opening, executing, and exporting the checked-in notebooks |

The old file also pinned unrelated application stacks such as Ansible, Dash,
Google authentication, and UMAP. No corresponding direct use was found in the
checked-in Python, notebooks, or shell commands. They are no longer direct
requirements; any packages needed transitively will still be installed by pip.

Local development now inherits the runtime's NumPy 2.3.4 and Matplotlib 3.11.2
pins, replacing the old environment capture's NumPy 2.1.3 and Matplotlib 3.10.0
pins. Keep an appropriate historical environment when reproducing older figures;
the curated developer environment follows the current pipeline runtime.

## Verification notes

The consolidated amd64 runtime passed all 62 regression tests, the pedigree
workflow with complete resume caching, and the trio smoke test. Developer
checks covered plotting, SVG-to-PDF conversion, Excel loading, notebook
execution and HTML export, and IGV-server/profiler imports in a native ARM64
Python 3.11 container.

Notebook kernel startup timed out in the amd64-emulated developer container on
macOS 26.6.2, while the same requirements passed natively. This is consistent
with a [reported Docker/Rosetta notebook issue](https://github.com/docker/desktop-feedback/issues/592).
Use the native local Python environment for notebook work; the pipeline CLI
checks succeeded under amd64 emulation. No Docker Desktop settings were changed
during verification.
