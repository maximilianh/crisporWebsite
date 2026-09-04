"""Create PRIDICT2.0 input sequences for flexible insertions/deletions.

This module contains the code that used to live in
``notebook_flexible_mutations.ipynb``.  The notebook imports from here, so
notebook and module cannot drift apart.

Typical module usage (no file I/O at all)::

    from flexible_mutation_input import flexible_mutation_sequences

    seqs = flexible_mutation_sequences(sequence, edit_type="insertion", insert_value="N")
    # -> [{'sequence_name': 'flexible_1', 'editseq': '...'}, ...]

The returned iterable can be handed straight to
``pridict2_pegRNA_design.predict_batch_sequences()`` or to
``parallel_batch_analysis(..., sequences=seqs)``.
"""

import itertools
import os

import pandas as pd

__all__ = [
    "IUPAC_CODES",
    "generate_edits",
    "flexible_mutation_sequences",
    "handle_duplicate_sequences",
    "flexible_mutation_input_generator",
]

# IUPAC nucleotide code mapping
IUPAC_CODES = {
    "A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"], "U": ["U"],
    "R": ["G", "A"], "Y": ["C", "T"], "K": ["G", "T"], "M": ["A", "C"],
    "S": ["G", "C"], "W": ["A", "T"], "B": ["G", "T", "C"], "D": ["G", "A", "T"],
    "H": ["A", "C", "T"], "V": ["G", "C", "A"], "N": ["A", "G", "C", "T"]
}


def generate_edits(sequence_name, sequence, edit_type, step=1, insert_value=None,
                   del_length=None, output_file_name=None, output_dir="./output",
                   visualise=False, verbose=True):
    """Generate all flexible insertion/deletion variants of ``sequence``.

    ``sequence`` has to contain the region in which the edit may be placed in
    square brackets, with at least 100 bp of context on both sides.

    Returns a dataframe with the columns ``sequence_name`` and ``editseq``.
    A csv file is only written when ``output_file_name`` is given.
    """
    # test if [] are in sequence, otherwise raise an error
    if "[" not in sequence or "]" not in sequence:
        raise ValueError("Sequence must contain brackets [] to indicate the target region of where to insert or delete.")
    # check if sequence_name, sequence, edit_type and step are provided and if either insert_value or del_length are provided
    if not sequence_name:
        raise ValueError("Sequence name not provided")
    if not sequence:
        raise ValueError("Sequence not provided")
    if edit_type not in ["insertion", "deletion"]:
        raise ValueError("Edit type must be either 'insertion' or 'deletion'")
    if not step:
        raise ValueError("Step size not provided")
    if step < 1:
        raise ValueError("Step size must be greater than 0")

    pre, target, post = sequence.split("[")[0], sequence.split("[")[1].split("]")[0], sequence.split("]")[1]
    edits = []

    # check if pre or post are < 100 bp long and if so, raise an error
    if len(pre) < 100 or len(post) < 100:
        error_parts = []
        if len(pre) < 100:
            error_parts.append(f"Sequence before '[' bracket is too short ({len(pre)} bp; minimum 100 bp required)")
        if len(post) < 100:
            error_parts.append(f"Sequence after ']' bracket is too short ({len(post)} bp; minimum 100 bp required)")
        raise ValueError(" and ".join(error_parts))

    if edit_type == "insertion":
        if len(target) < 1:
            raise ValueError("Target region is too short (0 bp; minimum 1 bp required to provide different options)")
        if insert_value is None:
            raise ValueError("Insertion sequence not provided.")
        # if any letter of edit_value is not a valid IUPAC code, return an empty list
        invalid_bases = [base for base in insert_value if base.upper() not in IUPAC_CODES]
        if invalid_bases:
            raise ValueError("Invalid DNA base(s): " + ", ".join(invalid_bases))
        possible_insertions = ["".join(p) for p in itertools.product(*[IUPAC_CODES.get(base.upper(), [base]) for base in insert_value])]
        for i in range(0, len(target) + 1, step):
            for insertion in possible_insertions:
                new_seq = pre + target[:i] + "(+" + insertion + ")" + target[i:] + post
                edits.append(new_seq)

    elif edit_type == "deletion":
        if len(target) < 2:
            raise ValueError("Target region is too short (1 bp; minimum 2 bp required to provide different deletion options.)")
        if del_length is None:
            raise ValueError("Deletion length not provided")
        elif del_length == 0:
            raise ValueError("Deletion length must be greater than 0 to be of any use")
        elif del_length >= len(target):
            raise ValueError("Deletion length must be shorter than the target region")
        if isinstance(del_length, float):  # if del_length is a float due to pandas import, convert to int if it is a whole number
            if del_length.is_integer():
                del_length = int(del_length)
            else:
                raise ValueError("Deletion length must be a whole integer")
        elif not isinstance(del_length, int):
            raise ValueError("Deletion length must be an integer")
        for i in range(0, len(target) - del_length + 1, step):
            new_seq = pre + target[:i] + "(-" + target[i:i+del_length] + ")" + target[i+del_length:] + post
            edits.append(new_seq)

    # create dataframe with columns sequence_name and editseq
    df = pd.DataFrame(columns=["sequence_name", "editseq"])
    df["sequence_name"] = [f"{sequence_name}_{i+1}" for i in range(len(edits))]
    df["editseq"] = edits

    # save dataframe to csv (only when an output filename is requested)
    if output_file_name:
        outpath = os.path.join(output_dir, output_file_name)
        os.makedirs(os.path.dirname(outpath) or ".", exist_ok=True)
        df.to_csv(outpath, index=False)
        if verbose:
            print(f"Saved {len(edits)} PRIDICT2 input sequences to {outpath}")

    if visualise:
        # Display the first 50 edits (only available inside a notebook)
        from IPython.display import display, HTML
        html_table = df.head(50).to_html()
        display(HTML("<div style='overflow-y: scroll; height: 300px;'>" + html_table + "</div>"))

    return df


def flexible_mutation_sequences(sequence, edit_type, step=1, insert_value=None,
                                del_length=None, sequence_name="flexible",
                                deduplicate=True, verbose=False):
    """Return an iterable of PRIDICT2.0 batch inputs for one target sequence.

    This is the module entry point: it takes a single input sequence (edit
    region in square brackets) and returns the rows that the notebook would
    have written to a csv file, as a list of
    ``{"sequence_name": ..., "editseq": ...}`` dictionaries.
    """
    df = generate_edits(
        sequence_name,
        sequence,
        edit_type,
        step=step,
        insert_value=insert_value,
        del_length=del_length,
        output_file_name=None,
        visualise=False,
        verbose=verbose,
    )
    if deduplicate:
        df = handle_duplicate_sequences(df, verbose=verbose)
    return df[["sequence_name", "editseq"]].to_dict("records")


def handle_duplicate_sequences(df, verbose=True):
    """
    Detect and handle duplicate sequence names in repetitive DNA sequences (identical edits), printing warnings and removing duplicates.
    Returns deduplicated dataframe.
    """
    # Find duplicate sequence names
    duplicates = df[df.duplicated(['sequence_name'], keep=False)]

    if not duplicates.empty:
        if verbose:
            # Group duplicates by sequence_name
            for name, group in duplicates.groupby('sequence_name'):
                print(f"\nIdentical sequence_name generated (resulting in the same edit):")
                print(f"Sequence name: {name}")
                print("PRIDICT inputs:")
                for idx, row in group.iterrows():
                    print(f"- {row['editseq']}")
                print("Removed duplicate sequences\n")

        # Keep first occurrence of each sequence_name
        df = df.drop_duplicates(subset=['sequence_name'], keep='first')

    return df


def flexible_mutation_input_generator(inputpath, inputfilename, outputpath=None,
                                      summaryoutputfilename=None,
                                      write_individual_files=False, verbose=True):
    """Batch mode: create flexible mutation inputs for every row of a csv file.

    Required columns: ``Name``, ``sequence_with_brackets``, ``edit_type``,
    ``insert``, ``deletion_length``, ``step``.

    Returns ``(input_df, final_df)``.  The summary csv is only written when
    both ``outputpath`` and ``summaryoutputfilename`` are given.
    """
    input_df = pd.read_csv(os.path.join(inputpath, inputfilename))
    all_dfs = []
    error_inputs = []
    for _, row in input_df.iterrows():
        sequence_name = row["Name"]
        sequence = row["sequence_with_brackets"]
        edit_type = row["edit_type"]
        insert = row["insert"]
        del_length = row["deletion_length"]
        step = row["step"]
        if write_individual_files:
            output_file_name = f"individual_{sequence_name}_flexible.csv"
        else:
            output_file_name = None

        if verbose:
            print(sequence_name)
        try:
            df = generate_edits(sequence_name, sequence, edit_type, step,
                                insert_value=insert, del_length=del_length,
                                output_file_name=output_file_name,
                                output_dir=outputpath or "./output",
                                visualise=False, verbose=verbose)
            all_dfs.append(df)
        except ValueError as e:
            if verbose:
                print(f"Error generating edits for {sequence_name}: {e}")
            error_inputs.append((sequence_name, e))
            continue

        if verbose:
            print()
    final_df = pd.concat(all_dfs, ignore_index=True)

    # Remove very rare duplicate bystander edit in repetitive sequence which are identical (e.g. CCT(CTC/tTt)TGG == C(CTC/tTt)TCTGG)
    final_df = handle_duplicate_sequences(final_df, verbose=verbose)

    if outputpath and summaryoutputfilename:
        summarypath = os.path.join(outputpath, summaryoutputfilename)
        os.makedirs(os.path.dirname(summarypath) or ".", exist_ok=True)
        final_df.to_csv(summarypath, index=False)
        if verbose:
            print(f"**********\nSaved {len(final_df)} PRIDICT2 input sequences to summary file {summarypath}.")

    if error_inputs and verbose:
        print("\nErrors occurred for the following sequences:")
        for seq_name, error in error_inputs:
            print(f"{seq_name}: {error}")
    return input_df, final_df
