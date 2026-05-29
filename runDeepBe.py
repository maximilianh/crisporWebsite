#!/data/www/crispor/venv/bin/python3
"""DeepBE: https://doi.org/10.1038/s41587-020-0573-5
https://github.com/MyungjaeSong/Paired-Library"""

import sys, os
import json
import importlib
from os.path import join, dirname

baseDir = join(dirname(__file__), "bin/deepBE")

modelList = json.load(open(join(baseDir, "modelList.json")))

# Work in progess

def run(data):
    """Runs the deepBE model. Returns dict with predicted efficiency
    and a list of (outcomeSeq, frequency). Placeholder for now."""

    # need to format data in crispor.py
    extGuideSeq = data[0]
    selModel = data[1]

    for model in modelList:
        if model == selModel:
            globals()[model] = importlib.import_module(model)

    return {"status": "processed", "output": data}
