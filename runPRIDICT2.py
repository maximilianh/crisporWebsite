"""Subserver entry point for pegRNA design with PRIDICT2.

The data posted by crispor.callSubServer("runPRIDICT2", data) is either

  * a plain string  - a PRIDICT2 input sequence, predicted as-is, or
  * a dict          - {"seq": <sequence>, "mode": <mode>, <mode options>}

"mode" selects how the sequence is turned into PRIDICT2 inputs:

  "single"           (default) predict the sequence exactly as given.
                     Sequence format: NNN(A/G)NNN, NNN(+TAG)NNN, NNN(-AC)NNN
                     with at least 100bp of context on both sides.

  "silentbystander"  design one variant per silent bystander edit around the
                     edit (bin/PRIDICT2/addons/silentbystander).
                     Sequence format: as for "single", but in frame and with
                     ~150bp of context on both sides - the bystander edits
                     widen the edit and therefore eat into the flanking context.
                     Options: silent ("yes"/"no"), changeEditBases ("no"/"yes"),
                     orfStart (0/1/2), surroundingAA (default 2)

  "flexibleedit"     design one variant per possible position of an insertion
                     or deletion (bin/PRIDICT2/addons/flexible_mutations).
                     Sequence format: the region in which the edit may be
                     placed in square brackets, e.g. NNN[ACGTACGT]NNN, with at
                     least 100bp of context on both sides of the brackets.
                     Options: editType ("insertion"/"deletion"), insert (the
                     inserted bases, IUPAC codes allowed), delLength, step

Options shared by all modes: name, maxSeqs, maxPegs, numProc.

Every mode returns the same structure, so a caller that only reads "out" does
not have to know which mode was used.
"""

import sys
from os.path import join, dirname

baseDir = dirname(__file__)
pridictDir = join(baseDir, "bin/PRIDICT2")
sys.path.insert(0, pridictDir)

from pridict2_pegRNA_design import predict_batch_sequences
from addons.flexible_mutations import flexible_mutation_sequences
from addons.silentbystander import silent_bystander_sequences

MODE_SINGLE = "single"
MODE_BYSTANDER = "silentbystander"
MODE_FLEXIBLE = "flexibleedit"

# accepted spellings of the "mode" option
MODES = {
    "": MODE_SINGLE,
    "none": MODE_SINGLE,
    "single": MODE_SINGLE,
    "standard": MODE_SINGLE,
    "silentbystander": MODE_BYSTANDER,
    "silent_bystander": MODE_BYSTANDER,
    "bystander": MODE_BYSTANDER,
    "flexibleedit": MODE_FLEXIBLE,
    "flexible_edit": MODE_FLEXIBLE,
    "flexible": MODE_FLEXIBLE,
    "flexiblemutation": MODE_FLEXIBLE,
    "flexible_mutations": MODE_FLEXIBLE,
}

# PRIDICT2 refuses sequences with less than 100bp flanking the edit
MIN_FLANK = 100

# Number of designed variants that are actually predicted. A variant costs ~9
# seconds, and crispor.py calls this subserver with a 60 second timeout, so the
# default is deliberately low: the addon modes design far more variants than
# that (a 1bp edit gives ~190 silent bystander variants). Raise "maxSeqs"
# together with the caller's timeout to score more of them.
DEFAULT_MAX_SEQS = 4
# One process is *faster* than several here: every extra process re-imports
# torch and re-loads the models, which costs more than the prediction itself.
DEFAULT_NUM_PROC = 1


def getOpt(data, *names, **kw):
    """Return the first of <names> present in data, else the default.

    Accepts camelCase and snake_case spellings of the same option."""
    default = kw.get("default")
    for name in names:
        if data.get(name) not in (None, ""):
            return data[name]
    return default


def flankLengths(editSeq):
    """Length of the sequence before and after the edit brackets.

    Handles both the PRIDICT format, NNN(A/G)NNN, and the flexibleedit input
    format, NNN[ACGT]NNN. Returns (-1, -1) if the sequence has no brackets."""
    for openBracket, closeBracket in (("(", ")"), ("[", "]")):
        if openBracket in editSeq and closeBracket in editSeq:
            return editSeq.index(openBracket), len(editSeq) - editSeq.index(closeBracket) - 1
    return -1, -1


def makeVariants(mode, seq, name, data):
    """Turn the input sequence into the list of PRIDICT2 inputs to predict.

    Returns a list of {"sequence_name", "editseq"} dicts; raises ValueError on
    invalid input."""
    if mode == MODE_SINGLE:
        if "(" not in seq or ")" not in seq:
            raise ValueError("Sequence must contain the edit in brackets, e.g. NNN(A/G)NNN.")
        return [{"sequence_name": name, "editseq": seq}]

    if mode == MODE_FLEXIBLE:
        editType = getOpt(data, "editType", "edit_type")
        if editType not in ("insertion", "deletion"):
            raise ValueError('"editType" must be either "insertion" or "deletion" in flexibleedit mode.')
        insert = getOpt(data, "insert", "insertValue", "insert_value")
        delLength = getOpt(data, "delLength", "deletionLength", "deletion_length")
        if delLength is not None:
            delLength = int(delLength)
        step = int(getOpt(data, "step", default=1))
        return flexible_mutation_sequences(
            seq,
            editType,
            step=step,
            insert_value=insert,
            del_length=delLength,
            sequence_name=name,
        )

    if mode == MODE_BYSTANDER:
        variants = silent_bystander_sequences(
            seq,
            name=name,
            silent=str(getOpt(data, "silent", default="yes")),
            change_edit_bases=str(getOpt(data, "changeEditBases", "change_edit_bases", default="no")),
            ORF_start=int(getOpt(data, "orfStart", "ORF_start", "orf_start", default=0)),
            silent_surrounding_AA_nr=int(
                getOpt(data, "surroundingAA", "silent_surrounding_AA_nr", default=2)
            ),
            full_records=True,
        )
        # when more variants are designed than can be predicted, keep the ones
        # that stay closest to the requested edit
        variants.sort(key=lambda v: (v["total_nr_of_base_changes"], v["sequence_name"]))
        return [{"sequence_name": v["sequence_name"], "editseq": v["editseq"]} for v in variants]

    raise ValueError("Unknown mode %r, expected one of: %s" % (mode, ", ".join(sorted(set(MODES.values())))))


def pegFromRow(pegDesc):
    """Build the peg description that crispor.showPegTable() expects."""
    oligos = [
        pegDesc["PCR-GG-Oligo1_Spacer"],
        pegDesc["PCR-GG-Oligo2_Extension"],
        pegDesc["Anza_Method_Spacer-Oligo-FW"],
        pegDesc["Anza_Method_Spacer-Oligo-RV"],
        pegDesc["Anza_Method_Extension-Oligo-FW"],
        pegDesc["Anza_Method_Extension-Oligo-RV"],
    ]
    return [
        pegDesc["pegRNA"],
        pegDesc["Spacer-Sequence"],
        pegDesc["PBSrevcomp"],
        pegDesc["RTrevcomp"],
        pegDesc["Target-Strand"],
        pegDesc["PRIDICT2_0_editing_Score_deep_K562"],
        pegDesc["PRIDICT2_0_editing_Score_deep_HEK"],
        pegDesc["Editing_Position"],
        pegDesc["protospacerlocation_only_initial"],
        pegDesc["PBSlocation"],
        pegDesc["RT_mutated_location"],
        pegDesc["Editor_Variant"],
        oligos,
    ]


def errorResult(mode, message):
    return {
        "status": "error",
        "model": "PRIDICT2",
        "mode": mode,
        "out": [],
        "variants": [],
        "log": message,
    }


def run(data):

    if isinstance(data, str):
        data = {"seq": data}
    if not isinstance(data, dict):
        return errorResult(MODE_SINGLE, "Expected a sequence string or a dict of options, got %s." % type(data).__name__)

    mode = str(getOpt(data, "mode", default=MODE_SINGLE)).strip().lower()
    if mode not in MODES:
        return errorResult(mode, "Unknown mode %r, expected one of: %s" % (mode, ", ".join(sorted(set(MODES.values())))))
    mode = MODES[mode]

    seq = getOpt(data, "seq", "editSeq", "sequence")
    if not seq:
        return errorResult(mode, 'No sequence given (expected the "seq" key).')
    seq = "".join(str(seq).split())

    name = str(getOpt(data, "name", "seqName", "sequence_name", default="test"))
    maxSeqs = int(getOpt(data, "maxSeqs", "max_seqs", default=DEFAULT_MAX_SEQS))
    maxPegs = getOpt(data, "maxPegs", "max_pegs")
    numProc = int(getOpt(data, "numProc", "num_proc", "cores", default=DEFAULT_NUM_PROC))

    try:
        variants = makeVariants(mode, seq, name, data)
    except ValueError as e:
        return errorResult(mode, str(e))
    except Exception as e:  # the addon modules also raise KeyError etc. on bad input
        return errorResult(mode, "Could not create %s input sequences: %s" % (mode, e))

    if not variants:
        return errorResult(mode, "No input sequences could be created for mode %s." % mode)

    # PRIDICT2 needs at least 100bp flanking the edit. The bystander edits widen
    # the edit, so variants of a sequence with just enough context can fall below.
    tooShort = [v for v in variants if min(flankLengths(v["editseq"])) < MIN_FLANK]
    variants = [v for v in variants if min(flankLengths(v["editseq"])) >= MIN_FLANK]
    if not variants:
        left, right = flankLengths(seq)
        return errorResult(
            mode,
            "Not enough flanking sequence: all %d designed variants have less than %dbp "
            "on one side of the edit (input has %dbp/%dbp). Provide more context."
            % (len(tooShort), MIN_FLANK, left, right),
        )

    nGenerated = len(variants) + len(tooShort)
    truncated = len(variants) > maxSeqs
    variants = variants[:maxSeqs]

    try:
        results = predict_batch_sequences(variants, num_proc=numProc)
    except Exception as e:
        return errorResult(mode, "pegRNA design failed: %s" % e)

    allPegs = []
    varInfos = []
    failed = []
    for variant in variants:
        varName = variant["sequence_name"]
        df = results.get(varName)
        if df is None:
            failed.append(varName)
            varInfos.append({"name": varName, "editseq": variant["editseq"], "pegCount": 0})
            continue
        pegs = [pegFromRow(pegDesc) for _, pegDesc in df.iterrows()]
        allPegs.extend(pegs)
        varInfos.append({"name": varName, "editseq": variant["editseq"], "pegCount": len(pegs)})

    if not allPegs:
        return errorResult(mode, "No pegRNA could be designed for any of the %d input sequences." % len(variants))

    # highest predicted K562 efficiency first, so that a caller which only keeps
    # the first pegRNAs keeps the best ones
    allPegs.sort(key=lambda peg: peg[5], reverse=True)
    if maxPegs:
        allPegs = allPegs[: int(maxPegs)]

    return {
        "status": "processed",
        "model": "PRIDICT2",
        "mode": mode,
        "out": allPegs,
        "variants": varInfos,
        "seqCount": {
            "generated": nGenerated,
            "skippedShort": len(tooShort),
            "predicted": len(variants),
            "failed": len(failed),
        },
        "truncated": truncated,
        "log": "",
    }
