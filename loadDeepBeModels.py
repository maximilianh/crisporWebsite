"""DeepBE: https://doi.org/10.1038/s41587-020-0573-5
https://github.com/MyungjaeSong/Paired-Library
This script loads all DeepBE models, which should get imported by runDeepBe.py later"""

import sys
import os
import json
import importlib
from os.path import join, dirname, isfile, getsize

baseDir = join(dirname(__file__), "bin/DeepBE")
jsonPath = join(baseDir, "modelList.json")

jsonData = []


def getModelDirs(modelDir):

    for subDir in os.listdir(modelDir):
        sys.path.append(join(modelDir, subDir))
        subPath = join(modelDir, subDir)
        if not os.path.isdir(subPath):
            continue
        for file in os.listdir(subPath):
            # get the name of modules
            if file.endswith(".py"):
                modelTpl = (subPath, file.strip(".py"))
                jsonData.append(modelTpl)


def loadAllModels():

    models = {}

    # load all models
    # models for SpCas9 fused with 7 deaminase domains
    SpCas9Models = join(baseDir, "DeepNG-BE")
    getModelDirs(SpCas9Models)

    # models for all combinations of 9 PAM variants / 7 deaminase domains
    pamVariantModels = join(baseDir, "DeepBE")
    getModelDirs(pamVariantModels)

    effModels = join(baseDir, "PAM")
    getModelDirs(effModels)

    if not isfile(jsonPath) or os.path.getsize(jsonPath) < 100:
        with open(jsonPath, "w", encoding="utf-8") as f:
            json.dump(jsonData, f)

    for _, name in json.load(open(jsonPath)):
        # load only some models for testing
        if "NG" not in name:
            continue

        mod = importlib.import_module(name)
        models[name] = mod.loadModel()   # captured + kept
        print("Imported %s" % name)

    return models
