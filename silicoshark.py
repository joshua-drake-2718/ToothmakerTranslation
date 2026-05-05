"""Top-level wrapper so `python silicoshark.py …` matches `python -m silicoshark …`.

Mirrors `main.py` (which runs Path A) for ergonomic parity. The actual
CLI lives in `silicoshark/__main__.py`; this file is a one-liner so
existing run-script invocations of the form

    python silicoshark.py PARAMS_FILE OUT_FOLDER OUT_NAME ITERATIONS SAVE_BLOCKS

work the same as the module form.
"""
from silicoshark.__main__ import main

if __name__ == '__main__':
    main()
