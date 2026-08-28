import os
import sys
from os.path import join, dirname

baseDir = dirname(__file__)
pridictDir = join(baseDir, "bin/PRIDICT2")
sys.path.insert(0, pridictDir)
outDir = join(pridictDir, "predictions")

from pridict2_pegRNA_design import predict_single_sequence


def run(data):

    seqName, editSeq = data

    editSeq = "GCCTGGAGGTGTCTGGGTCCCTCCCCCACCCGACTACTTCACTCTCTGTCCTCTCTGCCCAGGAGCCCAGGATGTGCGAGTTCAAGTGGCTACGGCCGA(G/C)GTGCGAGGCCAGCTCGGGGGCACCGTGGAGCTGCCGTGCCACCTGCTGCCACCTGTTCCTGGACTGTACATCTCCCTGGTGACCTGGCAGCGCCCAGATGCACCTGCGAACCACCAGAATGTGGCCGC"

    df = predict_single_sequence(
        sequence_name=seqName,
        editseq=editSeq,
        nicking=False,
        ngsprimer=False,
        use_5folds=False,
    )

    """
    allParams = "sequence_name,PRIDICT2_0_editing_Score_deep_K562,PRIDICT2_0_editing_Score_deep_HEK,K562_percentile_to_librarydiverse,HEK_percentile_to_librarydiverse,K562_rank,HEK_rank,PRIDICT2_Format,Original_Sequence,Edited_Sequence,Target-Strand,Mutation_Type,Correction_Type,Correction_Length,Editing_Position,PBSlength,RToverhanglength,RTlength,EditedAllele,OriginalAllele,Spacer-Sequence,PBSrevcomp,RTseqoverhangrevcomp,RTrevcomp,PCR-GG-Oligo1_Spacer,PCR-GG-Oligo2_Extension,Anza_Method_Spacer-Oligo-FW,Anza_Method_Spacer-Oligo-RV,Anza_Method_Extension-Oligo-FW,Anza_Method_Extension-Oligo-RV,Scaffold_Optimized,pegRNA,Editor_Variant,protospacermt,extensionmt,RTmt,RToverhangmt,PBSmt,original_base_mt,edited_base_mt,original_base_mt_nan,edited_base_mt_nan,RToverhangmatches,wide_initial_target,wide_mutated_target,protospacerlocation_only_initial,PBSlocation,RT_initial_location,RT_mutated_location,deepeditposition,deepeditposition_lst".split(",")
    """
    allPegs = []
    out = []
    maxRank = 0

    for _, pegDesc in df.iterrows():
        K562score = pegDesc[1]
        K562rank = pegDesc[3]
        HEKscore = pegDesc[2]
        # HEKPct = pegDesc[4]
        strand = pegDesc[10]
        editToNick = pegDesc[14]
        spacer = pegDesc[20]
        PBSrevComp = pegDesc[21]
        RTTrevComp = pegDesc[23]
        # primers
        oligos = [pegDesc[24], pegDesc[25],
                  pegDesc[26], pegDesc[27],
                  pegDesc[28], pegDesc[29]]

        pegSeq = pegDesc[31]
        editorVariant = pegDesc[32]
        spacerCoords = pegDesc[45]
        pbsCoords = pegDesc[46]
        rtCoords = pegDesc[48]

        allPegs.append((K562rank, [pegSeq, spacer, PBSrevComp, RTTrevComp, strand, K562score, HEKscore, editToNick, spacerCoords, pbsCoords, rtCoords, editorVariant, oligos]))
        maxRank = max(K562rank, maxRank)

    # only return the top 50 pegs with the highest K562 score
    for rank, peg in allPegs:
        """
        if rank < maxRank - 50:
            continue
        """
        out.append(peg)

    # placeholder
    return {"status": "processed",
            "model": "PRIDICT2",
            "out": out}
