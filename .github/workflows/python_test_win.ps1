python -m pip install pytest
if (-not $?) {throw "Failed to pip install pytest"}
python -m pip install fmm2dpy -f .\wheelhouse\
if (-not $?) {throw "Failed to pip install fmm2dpy"}
python -m pytest -s python/test
if (-not $?) {throw "Tests failed"}
