# PyCr2res

Python recipes for the [CRIRES+ pipeline](https://www.eso.org/sci/software/pipelines/cr2res/), built on ESO's pyesorex framework.

## Prerequisites

- **CPL library** (required): `apt-get install libcpl-dev` on Ubuntu/Debian
- **Python 3.12** (3.13 has compatibility issues)

## Quick start

```bash
git clone https://github.com/ivh/pycr2res
cd pycr2res
export PYESOREX_PLUGIN_DIR=$(pwd)/pyrecipes
uv sync
uv run pyesorex --recipes
```

## Structure

- `pyrecipes/` - Recipe files for pyesorex discovery
- `pycr2res/` - Python package with shared utilities

See `CLAUDE.md` for development details.
