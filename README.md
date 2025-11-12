# PyCr2res

This adds Python-based recipes to the [CRIRES+ pipeline](https://www.eso.org/sci/software/pipelines/cr2res/).

## Quick start

```bash
git clone https://github.com/ivh/pycr2res
cd pycr2res
export PYESOREX_PLUGIN_DIR=.
uv sync
uv run pyesorex --recipes
```

If you are building the whole pipeline, the appropriate place to put these files is `cr2rep/pyrecipes` but to get them installed together with the rest, you would first tell the build system about this directory.
