"""Create PRIDICT2.0 input sequences carrying silent (or any) bystander edits.

This module contains the code that used to live in
``notebook_silent_bystander_input.ipynb``.  The notebook imports from here, so
notebook and module cannot drift apart.

Typical module usage (no file I/O at all)::

    from silent_bystander_input import silent_bystander_sequences

    seqs = silent_bystander_sequences(pridict_input, name="cftr")
    # -> [{'sequence_name': 'cftr_G_c', 'editseq': '...'}, ...]

The returned iterable can be handed straight to
``pridict2_pegRNA_design.predict_batch_sequences()`` or to
``parallel_batch_analysis(..., sequences=seqs)``.
"""

import math
import os
import re
from itertools import product

import pandas as pd
from Bio.Seq import Seq

__all__ = [
    "generate_all_sequences",
    "generate_combinations",
    "convert_differences_to_lowercase",
    "split_sequence",
    "validate_context_length",
    "process_contexts",
    "handle_duplicate_sequences",
    "isDNA",
    "primesequenceparsing",
    "bystander_creation_for_pridict",
    "silent_bystander_sequences",
    "bystander_input_generator",
]

# default variables (do not change unless you want to adapt the script)
# number of amino acids up- and downstream of edit AA for which silent bystander
# mutations will be created; default = 2; maximum 4 AA are allowed in this script
DEFAULT_SILENT_SURROUNDING_AA_NR = 2
# limit maximum length of edit (including bystander edits) to 40 bases;
# max. edit length for PRIDICT2.0 predictions is 40 bases
DEFAULT_TOTAL_EDIT_LIMIT = 40
# maximum length of the edit bases you want to change
DEFAULT_MAX_EDIT_LENGTH = 10
# minimum edit-flanking length, after correcting for ORF_start. Only required as sanity check.
DEFAULT_MINIMUM_FLANKING = 94


def generate_all_sequences(length):
    ### Generate all possible sequences of a given length
    nucleotides = ['A', 'T', 'G', 'C']
    return [''.join(seq) for seq in product(nucleotides, repeat=length)]

def generate_combinations(left_context_close, left_options, edit, right_options, right_context_close):
    ### Generate all possible combinations of left and right context with edit
    return [
        left_context_close + left_opt + edit + right_opt + right_context_close
        for left_opt, right_opt in product(left_options, right_options)
    ]

def convert_differences_to_lowercase(option, original):
    ### Convert differences between option and original to lowercase
    return ''.join(
        opt.lower() if opt != orig else opt 
        for opt, orig in zip(option, original)
    )

def split_sequence(seq):
    ### Split a sequence into three parts: before first lowercase letter, between first and last lowercase letter, and after last lowercase letter
    lowercase_positions = [m.start() for m in re.finditer('[a-z]', seq)]
    if not lowercase_positions:
        return seq, '', '', 0
    
    first_lower_pos = lowercase_positions[0]
    last_lower_pos = lowercase_positions[-1]
    before_first_lower = seq[:first_lower_pos]
    between_first_last_lower = seq[first_lower_pos:last_lower_pos + 1]
    after_last_lower = seq[last_lower_pos + 1:]

    return before_first_lower, between_first_last_lower, after_last_lower, len(lowercase_positions)

def validate_context_length(input, minimum_flanking):
    ### Validate that the context length is at least the minimum flanking length
    left_context = len(input.split('(')[0])
    right_context = len(input.split(')')[1])
    if left_context < minimum_flanking or right_context < minimum_flanking:
        raise ValueError('Context length is less than minimum flanking length! Please check your input sequence.')

def process_contexts(input, bystander_window, close_context_len):
    ### Process the input sequence into left and right context
    left_context = input.split('(')[0]
    right_context = input.split(')')[1]
    left_context_start = left_context[:-bystander_window-close_context_len]
    left_context_close = left_context[-bystander_window-close_context_len:-bystander_window]
    left_bystander_window = left_context[-bystander_window:]
    right_context_end = right_context[bystander_window+close_context_len:]
    right_context_close = right_context[bystander_window:bystander_window+close_context_len]
    right_bystander_window = right_context[:bystander_window]
    return left_context_start, left_context_close, left_bystander_window, right_context_close, right_bystander_window, right_context_end

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
                print("Removed additional bystander variants\n")
        
        # Keep first occurrence of each sequence_name
        df = df.drop_duplicates(subset=['sequence_name'], keep='first')
    
    return df

def isDNA(sequence):
    """ Check whether sequence contains only DNA bases. """
    onlyDNA = True
    diff_set = set(sequence) - set('ACTGatgc')
    if diff_set:
        onlyDNA = False
        print('Non-DNA bases detected. Please use ATGC.')
        print(sequence)
        raise ValueError
    return onlyDNA

def primesequenceparsing(sequence: str) -> object:
    """
    Function which takes target sequence with desired edit as input and 
    editing characteristics as output. Edit within brackets () and original
    equence separated by backslash from edited sequence: (A/G) == A to G mutation.
    Placeholder for deletions and insertions is '-'.

    Parameters
    ----------
    sequence : str
        Target sequence with desired edit in brackets ().
    """
    
    sequence = sequence.replace('\n','')  # remove any spaces or linebreaks in input
    sequence = sequence.replace(' ','')
    sequence = sequence.upper()
    if sequence.count('(') != 1:
        print(sequence)
        print('More or less than one bracket found in sequence! Please check your input sequence.')
        raise ValueError

    five_prime_seq = sequence.split('(')[0]
    three_prime_seq = sequence.split(')')[1]

    sequence_set = set(sequence)
    if '/' in sequence_set:
        original_base = sequence.split('/')[0].split('(')[1]
        edited_base = sequence.split('/')[1].split(')')[0]

        # edit flanking bases should *not* be included in the brackets
        if (original_base[0] == edited_base[0]) or (original_base[-1] == edited_base[-1]):
            print(sequence)
            print('Flanking bases should not be included in brackets! Please check your input sequence.')
            raise ValueError
    elif '+' in sequence_set:  #insertion
        original_base = '-'
        edited_base = sequence.split('+')[1].split(')')[0]
    elif '-' in sequence_set:  #deletion
        original_base = sequence.split('-')[1].split(')')[0]
        edited_base = '-'

    # ignore "-" in final sequences (deletions or insertions)
    if original_base == '-':
        original_seq = five_prime_seq + three_prime_seq
        if edited_base != '-':
            mutation_type = 'Insertion'
            correction_length = len(edited_base)
        else:
            print(sequence)
            raise ValueError
    else:
        original_seq = five_prime_seq + original_base + three_prime_seq
        if edited_base == '-':
            mutation_type = 'Deletion'
            correction_length = len(original_base)
        elif len(original_base) == 1 and len(edited_base) == 1:
            if isDNA(original_base) and isDNA(edited_base):  # check if only AGCT is in bases
                mutation_type = '1bpReplacement'
                correction_length = len(original_base)
            else:
                print(sequence)
                print('Non DNA bases found in sequence! Please check your input sequence.')
                raise ValueError
        elif len(original_base) > 1 or len(edited_base) > 1:
            if isDNA(original_base) and isDNA(edited_base):  # check if only AGCT is in bases
                mutation_type = 'MultibpReplacement'
                if len(original_base) == len(
                        edited_base):  # only calculate correction length if replacement does not contain insertion/deletion
                    correction_length = len(original_base)
                else:
                    print(sequence)
                    print('Only 1bp replacements or replacements of equal length (before edit/after edit) are supported! Please check your input sequence.')
                    raise ValueError
            else:
                print(sequence)
                print('Non DNA bases found in sequence! Please check your input sequence.')
                raise ValueError

    if edited_base == '-':
        edited_seq = five_prime_seq + three_prime_seq
    else:
        edited_seq = five_prime_seq + edited_base.lower() + three_prime_seq

    if isDNA(edited_seq) and isDNA(original_seq):  # check whether sequences only contain AGCT
        pass
    else:
        raise ValueError

    basebefore_temp = five_prime_seq[
                      -1:]  # base before the edit, could be changed with baseafter_temp if Rv strand is targeted (therefore the "temp" attribute)
    baseafter_temp = three_prime_seq[:1]  # base after the edit

    editposition_left = len(five_prime_seq)
    editposition_right = len(three_prime_seq)
    return original_base, edited_base, original_seq, edited_seq, editposition_left, editposition_right, mutation_type, correction_length, basebefore_temp, baseafter_temp

def bystander_creation_for_pridict(pridict_input_original, silent_surrounding_AA_nr, ORF_start, name,
                                   minimum_flanking, total_edit_limit, max_edit_length,
                                   silent='yes', change_edit_bases='no', verbose=True):
    if silent == 'no':
        if silent_surrounding_AA_nr != 1:
            if verbose:
                print('*** Flanking AA number is set to 1 (maximum) for non-silent bystander, as more will lead to memory issues due to a huge amount of possibilities. ***')
            silent_surrounding_AA_nr = 1



    # Map each amino acid to its corresponding codons
    codon_map = {
        'A': ['GCT', 'GCC', 'GCA', 'GCG'],
        'C': ['TGT', 'TGC'],
        'D': ['GAT', 'GAC'],
        'E': ['GAA', 'GAG'],
        'F': ['TTT', 'TTC'],
        'G': ['GGT', 'GGC', 'GGA', 'GGG'],
        'H': ['CAT', 'CAC'],
        'I': ['ATT', 'ATC', 'ATA'],
        'K': ['AAA', 'AAG'],
        'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
        'M': ['ATG'],
        'N': ['AAT', 'AAC'],
        'P': ['CCT', 'CCC', 'CCA', 'CCG'],
        'Q': ['CAA', 'CAG'],
        'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
        'T': ['ACT', 'ACC', 'ACA', 'ACG'],
        'V': ['GTT', 'GTC', 'GTA', 'GTG'],
        'W': ['TGG'],
        'Y': ['TAT', 'TAC'],
        '*': ['TAA', 'TAG', 'TGA']  # Stop codons
    }
    # assert that ORF_start is in [0, 1, 2]
    if ORF_start not in [0, 1, 2]:
        raise ValueError('ORF_start has to be 0, 1 or 2!')
    
    pridict_input_original = pridict_input_original[ORF_start:]
    minimum_flanking = minimum_flanking - ORF_start # subtract ORF_start from minimum_flanking to account for the shift in the ORF and allow that the minimum flanking is still met
    pridict_input_original = pridict_input_original.upper()
    validate_context_length(pridict_input_original, minimum_flanking)

    # confirm that pridict_input_original contains replacement edit:
    mutation_type = primesequenceparsing(pridict_input_original)[6]
    if not mutation_type in ['1bpReplacement', 'MultibpReplacement']:
        raise ValueError('Only 1bp replacements or multi bp replacements of equal length (before edit/after edit) are supported for silent bystander predictions! Please check your input sequence.')

    left_context = pridict_input_original.split('(')[0]  # bases before brackets
    right_context = pridict_input_original.split(')')[1]  # bases after brackets
    left_context_length = len(left_context)
    # number of bases which are in the context, but already part of the AA that is edited
    remaining_bases = left_context_length % 3

    # sequence/base of the defined target before editing:
    edit_before = pridict_input_original.split('(')[1].split('/')[0]
    # sequence/base of the defined target after editing:
    edit_after = pridict_input_original.split('/')[1].split(')')[0]
    edit_after = convert_differences_to_lowercase(edit_after, edit_before)  # converts the actual edited bases to lowercase; not the whole edited sequence (e.g. AGA to TGT will be tGt and not tgt)
    edit_length = len(edit_before)
    if edit_length > max_edit_length:
        raise ValueError(f'Only edits of up to {max_edit_length} bases are supported! Please check your input sequence.')

    # sequence of the edited target including context
    edited_sequence = left_context + edit_after + right_context
    # sequence of the original target including context
    original_sequence = left_context + edit_before + right_context

    # number of bases which make up the AA(s) that are part of the edit (in brackets of the pridict input)
    edit_AA_length_in_bp = math.ceil((remaining_bases + edit_length)/3)*3
    # number of AA that are edited (silent or not)
    edit_AA_length = int(edit_AA_length_in_bp / 3)
    downstream_AA_context_basenr = edit_AA_length_in_bp - edit_length - remaining_bases

    # sequence upstream that is unchanged (also not silently changed)
    untouched_left_context = left_context[:-remaining_bases-silent_surrounding_AA_nr*3]
    # sequence downstream that is unchanged (also not silently changed)
    untouched_right_context = right_context[downstream_AA_context_basenr+silent_surrounding_AA_nr*3:]

    # sequence of bases that can be changed by the silent bystander mutations (including edit bases)
    start_pos = len(left_context)-remaining_bases-silent_surrounding_AA_nr*3
    end_pos = len(left_context)-remaining_bases+edit_AA_length_in_bp+silent_surrounding_AA_nr*3

    original_bases = original_sequence[start_pos:end_pos]
    original_AA_seq = str(Seq(original_bases).translate())

    potentially_changed_bases = edited_sequence[start_pos:end_pos]
    potentially_changed_AA = str(Seq(potentially_changed_bases).translate())

    if silent == 'yes':
    # create all possible codon variants which result in the same AA sequence
        bystander_option_list = []
        for aa in potentially_changed_AA:
            amino_acid_codons = codon_map[aa]
            bystander_option_list.append(amino_acid_codons)

        # generate all possible combinations of codons in the bystander_option_list. For each AA, we have a list of possible codons in there
        all_combinations = list(product(*bystander_option_list))
        result_strings = [''.join(combination) for combination in all_combinations]

        # remove the mutated sequence without any bystander from this list
        result_strings = [result_string for result_string in result_strings if result_string.upper() != potentially_changed_bases.upper()]
        if verbose:
            print('Number of silent options INCLUDING changed edit bases', len(result_strings))

        if change_edit_bases == 'no':
            # remove options where also defined edited bases are changed (silently); check whether every position where we have a lowercase letter in potentially_changed_bases. this position should be the same base in the result_string and potentially_changed_bases
            result_strings = [result_string for result_string in result_strings if all([result_string[i].upper() == potentially_changed_bases[i].upper() for i in range(len(result_string)) if potentially_changed_bases[i].islower()])]
            if verbose:
                print('Number of silent options EXCLUDING changed edit bases', len(result_strings))

        result_AAs = [str(Seq(result_string).translate()) for result_string in result_strings]


    elif silent == 'no':
        result_strings = generate_all_sequences(len(potentially_changed_bases))

        # remove the mutated sequence without any bystander from this list
        result_strings = [result_string for result_string in result_strings if result_string.upper() != potentially_changed_bases.upper()]
        if verbose:
            print('Number of bystander options INCLUDING changed edit bases', len(result_strings))

        if change_edit_bases == 'no':
            # remove options where also defined edited bases are changed (silently); check whether every position where we have a lowercase letter in potentially_changed_bases. this position should be the same base in the result_string and potentially_changed_bases
            result_strings = [result_string for result_string in result_strings if all([result_string[i].upper() == potentially_changed_bases[i].upper() for i in range(len(result_string)) if potentially_changed_bases[i].islower()])]
            if verbose:
                print('Number of bystander options EXCLUDING changed edit bases', len(result_strings))

        result_AAs = [str(Seq(result_string).translate()) for result_string in result_strings]

    # for each element of result_string, check whether a base is different from the original sequence. If so, convert it to lowercase
    result_strings = [convert_differences_to_lowercase(result_string, original_bases) for result_string in result_strings]

    full_edited_list = [untouched_left_context + result_string + untouched_right_context for result_string in result_strings]

    pridict_input_lst = []
    totalbasechanges_lst = []
    edit_name_bystander_lst = []
    original_length_lst = []
    final_length_lst = []
    edited_bystander_focus_sequence_lst = []
    edited_bystander_focus_AA_lst = []
    editonly_focus_sequence_lst = []
    editonly_nobystander_focus_AA_lst = []
    for index, sequence in enumerate(full_edited_list):
        unchanged_before, edit, unchanged_after, totalbasechanges = split_sequence(sequence)
        original_base = original_sequence[len(unchanged_before):-len(unchanged_after)]
        final_pridict_seq = f"{unchanged_before}({original_base}/{edit}){unchanged_after}"
        edit_name_bystander = f"{name}_{original_base}_{edit}"
        edited_bystander_focuse_sequence = result_strings[index]
        edited_bystander_focus_AA = result_AAs[index]
        edit_name_bystander_lst.append(edit_name_bystander)
        final_length_lst.append(len(edit))
        pridict_input_lst.append(final_pridict_seq)
        totalbasechanges_lst.append(totalbasechanges)
        original_length_lst.append(edit_length)
        edited_bystander_focus_sequence_lst.append(edited_bystander_focuse_sequence)
        edited_bystander_focus_AA_lst.append(edited_bystander_focus_AA)
        editonly_focus_sequence_lst.append(potentially_changed_bases)
        editonly_nobystander_focus_AA_lst.append(potentially_changed_AA)

    if verbose:
        print(f'Number of possible bystander mutations: {len(result_strings)}')
        print('Sequence flanking edit position (before edit):',original_bases)
        print('Sequence flanking edit position (AFTER edit, without bystander editing):',potentially_changed_bases)
        print('AA flanking edit position (before edit):',original_AA_seq)
        print('AA flanking edit position (after edit, without bystander editing):',potentially_changed_AA)
        print('AA flanking edit position (after edit, WITH bystander editing:',set(result_AAs))

    finaldf = pd.DataFrame(list(zip(
        edit_name_bystander_lst, pridict_input_lst, original_length_lst, final_length_lst, totalbasechanges_lst, edited_bystander_focus_sequence_lst, edited_bystander_focus_AA_lst, editonly_focus_sequence_lst, editonly_nobystander_focus_AA_lst
    )), columns=['sequence_name', 'editseq', 'original_edit_length', 'final_edit_length_with_bystander', 'total_nr_of_base_changes', 'bystander_focus_sequence', 'bystander_focus_AA', 'editedonly_focus_sequence', 'editedonly_focus_AA'])

    # filter finaldf to exclude rows where final_edit_length_with_bystander > total_edit_limit
    finaldf = finaldf[finaldf.final_edit_length_with_bystander <= total_edit_limit]
    return finaldf


def silent_bystander_sequences(pridict_input, name="bystander", silent="yes",
                               change_edit_bases="no", ORF_start=0,
                               silent_surrounding_AA_nr=DEFAULT_SILENT_SURROUNDING_AA_NR,
                               total_edit_limit=DEFAULT_TOTAL_EDIT_LIMIT,
                               max_edit_length=DEFAULT_MAX_EDIT_LENGTH,
                               minimum_flanking=DEFAULT_MINIMUM_FLANKING,
                               deduplicate=True, full_records=False, verbose=False):
    """Return an iterable of PRIDICT2.0 batch inputs for one target sequence.

    This is the module entry point: it takes a single PRIDICT input sequence
    (with 150 bp context on both sides of the edit) and returns the rows that
    the notebook would have written to a csv file, as a list of dictionaries.

    With ``full_records=False`` (default) each dict holds only
    ``sequence_name`` and ``editseq`` - the two columns PRIDICT2.0 batch mode
    needs.  Set ``full_records=True`` to also get the bystander annotation
    columns.
    """
    df = bystander_creation_for_pridict(
        pridict_input,
        silent_surrounding_AA_nr,
        ORF_start,
        name,
        minimum_flanking,
        total_edit_limit,
        max_edit_length,
        silent=silent,
        change_edit_bases=change_edit_bases,
        verbose=verbose,
    )
    if deduplicate:
        df = handle_duplicate_sequences(df, verbose=verbose)
    if not full_records:
        df = df[["sequence_name", "editseq"]]
    return df.to_dict("records")


def bystander_input_generator(inputpath, inputfilename, outputpath=None, outputfilename=None,
                              silent_surrounding_AA_nr=DEFAULT_SILENT_SURROUNDING_AA_NR,
                              total_edit_limit=DEFAULT_TOTAL_EDIT_LIMIT,
                              max_edit_length=DEFAULT_MAX_EDIT_LENGTH,
                              minimum_flanking=DEFAULT_MINIMUM_FLANKING, verbose=True):
    """Batch mode: create bystander inputs for every row of a csv file.

    Required columns: ``Name``, ``pridict_input``, ``silent``,
    ``change_edit_bases``, ``in_frame``.

    Returns ``(input_df, final_df)``.  The output csv is only written when both
    ``outputpath`` and ``outputfilename`` are given.
    """
    df = pd.read_csv(os.path.join(inputpath, inputfilename))
    all_dfs = []
    for _, row in df.iterrows():
        name = row["Name"]
        silent = row["silent"]
        change_edit_bases = row["change_edit_bases"]
        if verbose:
            print(name)
        ORF_start = row["in_frame"]
        if ORF_start == 'yes':
            ORF_start = 0
        else:
            raise ValueError('ORF start is not in frame! Please check your input sequence. Input sequence has to be in-frame.')
        pridict_input_original = row.pridict_input
        silentbystanderdf = bystander_creation_for_pridict(pridict_input_original, silent_surrounding_AA_nr, ORF_start, name, minimum_flanking, total_edit_limit, max_edit_length, silent=silent, change_edit_bases=change_edit_bases, verbose=verbose)
        if verbose:
            print(f'{len(silentbystanderdf)} silent bystander sequences created for {name}')
            print()
        all_dfs.append(silentbystanderdf)
    final_df = pd.concat(all_dfs, ignore_index=True)

    # Remove very rare duplicate bystander edit in repetitive sequence which are identical (e.g. CCT(CTC/tTt)TGG == C(CTC/tTt)TCTGG)
    final_df = handle_duplicate_sequences(final_df, verbose=verbose)

    if outputpath and outputfilename:
        outfile = os.path.join(outputpath, outputfilename)
        os.makedirs(os.path.dirname(outfile) or ".", exist_ok=True)
        final_df.to_csv(outfile)
        if verbose:
            print(f"Saved {len(final_df)} PRIDICT2 input sequences to {outfile}.")
    return df, final_df
