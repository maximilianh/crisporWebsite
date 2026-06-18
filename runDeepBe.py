"""DeepBE: https://doi.org/10.1038/s41587-020-0573-5
https://github.com/MyungjaeSong/Paired-Library"""

import loadDeepBeModels
import importlib
from os.path import join, dirname

baseDir = join(dirname(__file__), "bin/deepBE")
MODELS = loadDeepBeModels.loadAllModels()   # runs ONCE, at subserver startup

pamToModel = {
        0: 'PAM_variant_SpCas9',
        1: 'PAM_variant_VRQR',
        2: 'PAM_variant_NG',
        3: 'PAM_variant_NRRH',
        4: 'PAM_variant_NRTH',
        5: 'PAM_variant_NRCH',
        6: 'PAM_variant_SpG',
        7: 'PAM_variant_SpRY',
        8: 'PAM_variant_Sc++'
        }


def run(data):
    """Runs the deepBE model. Returns dict with predicted efficiency
    and a list of (outcomeSeq, frequency). Placeholder for now."""

    # Work in progess
    editor, selModel, pamVariant, extGuideSeq = data

    outcomeModule = importlib.import_module(selModel)

    if selModel[0:5] == "DeepNG":
        pamVariant = 2
        mainModel = MODELS[selModel]
        outcomes = outcomeModule.predict(mainModel, [extGuideSeq], pamVariant)

    else:
        # DeepBE needs models for PAM variants too
        # need to select the PAM variant (info stored in DeepBE/PamModelId.json
        # temporary placeholder
        pamVariant = 2
        mainModel = MODELS[selModel]
        outcomes = outcomeModule.predict(mainModel, [extGuideSeq], pamVariant)

    # import Cas9 efficiency model
    """
    effModule = importlib.import_module(pamToModel[pamVariant])
    effModel = MODELS[pamToModel[pamVariant]]
    effs = effModule.predict(effModel, [extGuideSeq], pamVariant)

    for _, effRow in effs.iterrows():
        # divide by 100 to match the processing of other scores in crispor.py
        totalEdit = float(effRow["Prediction score"])
        break

    """

    # get efficiency as the sum of outcome frequencies (total editing)
    totalEdit = 0

    # transform the pandas data frame into a list
    outcomeList = []
    for _, outcomeRow in outcomes.iterrows():
        # make edited bases uppercase and the rest in lowercase
        freq = float(outcomeRow["Prediction score"])
        outcomeList.append((outcomeRow["edited output"], freq))
        totalEdit += freq

    return {"status": "processed",
            "model": selModel,
            "eff": totalEdit,
            "outcome": outcomeList}
