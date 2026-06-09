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
    editor, selModel, pamVariant, extGuideSeq = data

    mod = importlib.import_module(selModel)

    if selModel[0:5] == "DeepNG":
        pamVariant = 2
        mainModel = MODELS[selModel]
        outcomes = mod.predict(mainModel, [extGuideSeq], pamVariant)   # reuse in-memory model

    else:
        # DeepBE needs models for PAM variants too
        # need to select the PAM variant (info stored in DeepBE/PamModelId.json
        # temporary placeholder
        pamVariant = 2
        mainModel = MODELS[selModel]
        outcomes = mod.predict(mainModel, [extGuideSeq], pamVariant)  # reuse in-memory model

    # get efficiency as the sum of outcome frequencies (total editing)
    totalEff = 0

    # transform the pandas data frame into a list
    outcomeList = []
    for _, row in outcomes.iterrows():
        # make edited bases uppercase and the rest in lowercase
        freq = float(row["Prediction score"])
        outcomeList.append((row["edited output"], freq))
        totalEff += freq

    return {"status": "processed",
            "model": selModel,
            "eff": totalEff,
            "outcome": outcomeList}
