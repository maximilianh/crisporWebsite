import sys

sys.path.append("bin/FORECasT-BE")
import forecast_be as forecast


def run(data):

    """
    calcultes the predicted efficiency and outcomes with FORECast-BE
    see https://github.com/ananth-pallaseni/FORECasT-BE
    """


    # forecast.load_models()

    extGuideSeq = data
    editor = "CBE"
    # do editor = data[1]

    # to match the output of deepBE / CRISPRonBE, the 30bp extended guide sequence is returned
    guideSeq = extGuideSeq[4:24]

    # Predict the total fraction of edited reads for the target sequence
    # If the `mean` and `std` arguments are None, then returns a z-score
    # Input a mean and std to scale this into reael efficiency (good defaults are mean=0.5 & std=0.1)
    mean, std = 0.5, 0.1

    totalEff = forecast.predict_total(guideSeq, editor=editor, mean=mean, std=std)

    # Predict the fraction of edited reads with the on-target substitituion at each position
    # Returns a list of predictions
    outcomes = []
    posEff = forecast.predict(guideSeq, editor=editor, mean=mean, std=std)
    for pos, freq in posEff:
        if freq is None:
            continue
        idx = pos - 1
        # will need to adjust for ABE and reverse strand using fromNucl / toNucl
        outcomeSeq = extGuideSeq[0:4].lower() + guideSeq[0:idx].lower() + "T" + guideSeq[pos:].lower() + extGuideSeq[24:].lower()
        outcomes.append((outcomeSeq, freq))

        # print(outcomes, freq, "<br>")

    return {"status": "processed",
            "model": "FORECast-BE",
            "eff": totalEff,
            "outcome": outcomes}
