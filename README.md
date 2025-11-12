# PyCr2res

This adds Python-based recipes to the [CRIRES+ pipeline](https://www.eso.org/sci/software/pipelines/cr2res/).

Simply set PYESOREX_PLUGIN_DIR to the directory with these files, the `pyesorex --recipes` should show the new recipes.

If you are building the whole pipeline, the appropriate place to put the files is `cr2rep/pyrecipes` but to get the files installed together with the rest, you would first tell the build system about them.
