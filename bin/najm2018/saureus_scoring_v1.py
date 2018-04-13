import pandas as pd
import csv, argparse
import azimuth
import azimuth.model_comparison
import numpy as np

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-file',
        type=str,
        help='Input file with context')
    parser.add_argument('--model-file',
        type=str,
        help='Model file')
    parser.add_argument('--outputfile',
        type=str,
        help='Outputfile')
    return parser

if __name__ == '__main__':
    args = get_parser().parse_args()
    input_file = args.input_file
    outputfile = args.outputfile
    model_file = args.model_file
    input_df = pd.read_table(input_file)
    colnames = list(input_df.columns)
    colnames.append('On-Target Efficacy Score')
    with open(outputfile,'w') as o:
        w = csv.writer(o,delimiter='\t',lineterminator='\n')
        w.writerow(colnames)
        for i,r in input_df.iterrows():
            if len(r[1]) == 36:
                score = azimuth.model_comparison.predict(np.array([r[1].upper()]), aa_cut=None, percent_peptide=None, model_file=model_file, pam_audit=False)
                row = list(r)
                row.append(score[0])
                w.writerow(row)
            else:
                sys.exit("Context sequence should be 36 nucleotides long")

