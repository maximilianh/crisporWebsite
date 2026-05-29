#!/data/www/crispor/venv/bin/python3
"""DeepBE: https://doi.org/10.1038/s41587-020-0573-5
https://github.com/MyungjaeSong/Paired-Library
This script loads all DeepBE models, which should get imported by runDeepBe.py later"""

import sys
import os
import json
import importlib
from os.path import join, dirname

baseDir = join(dirname(__file__), "bin/DeepBE")

# optionally write a JSON file to store the name and path of all modules instead of loading them
writeJson = False
jsonData = []


def getModelDirs(modelDir, writeJson):

    for subDir in os.listdir(modelDir):
        sys.path.append(join(modelDir, subDir))
        subPath = join(modelDir, subDir)
        if not os.path.isdir(subPath) or not writeJson:
            continue
        for file in os.listdir(subPath):
            # get the name of modules
            if file.endswith(".py"):
                modelTpl = (subPath, file.strip(".py"))
                jsonData.append(modelTpl)


def main():

    # load all models
    # models for SpCas9 fused with 7 deaminase domains
    SpCas9Models = join(baseDir, "DeepNG-BE")
    getModelDirs(SpCas9Models, writeJson)

    # models for all combinations of 9 PAM variants / 7 deaminase domains
    pamVariantModels = join(baseDir, "DeepBE")
    getModelDirs(pamVariantModels, writeJson)

    jsonPath = join(baseDir, "modelList.json")

    if writeJson:
        with open(jsonPath, "w", encoding="utf-8") as f:
            json.dump(jsonData, f)
    else:
        modelList = json.load(open(jsonPath))
        for modelTpl in modelList:
            model = modelTpl[1]
            globals()[model] = importlib.import_module(model)


if __name__ == "__main__":
    main()
