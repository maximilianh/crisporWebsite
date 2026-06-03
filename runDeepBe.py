"""DeepBE: https://doi.org/10.1038/s41587-020-0573-5
https://github.com/MyungjaeSong/Paired-Library"""

import loadDeepBeModels
import importlib
from os.path import join, dirname

baseDir = join(dirname(__file__), "bin/deepBE")
MODELS = loadDeepBeModels.loadAllModels()   # runs ONCE, at subserver startup


def run(data):
    """Runs the deepBE model. Returns dict with predicted efficiency
    and a list of (outcomeSeq, frequency). Placeholder for now."""

    # Work in progess
    editor, extGuideSeq = data

    if editor == "CBE":
        # select the model
        selModel = "DeepBE_CGBE1"

    mod = importlib.import_module(selModel)
    # no effs ?
    totalEff = 0

    if selModel[0:5] == "DeepNG":
        pamVariant = 2
        mainModel = MODELS[selModel]
        outcomes = mod.predict(mainModel, list(extGuideSeq), pamVariant)   # reuse in-memory model

    else:
        # DeepBE needs models for PAM variants too
        # need to select the PAM variant (info stored in DeepBE/PamModelId.json
        # temporary placeholder
        pamVariant = 2
        mainModel, pamModels = MODELS[selModel]
        print(mainModel)
        outcomes = mod.predict(mainModel, pamModels, list(extGuideSeq), pamVariant)   # reuse in-memory model
        print(outcomes)

        pass

    outcomes = ["AAA"]

    return {"status": "processed",
            "model": selModel,
            "eff": totalEff,
            "outcome": outcomes}
