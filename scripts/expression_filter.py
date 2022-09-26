import os
import argparse
import pandas as pd

def parse_args():
    description = "Filter neoantigen lines with Kallisto TPM"
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--pvacseq", help="pVACseq neoantigen filtered.tsv")
    p.add_argument("--kallisto", help="Kallisto abundance.tsv")
    p.add_argument("--cutoff", help="Minimum TPM cutoff")
    p.add_argument("--output", help="pVACseq neoantigen filtered by TPM")
    args = p.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    assert os.path.exists(args.pvacseq), args.pvacseq
    assert os.path.exists(args.kallisto), args.kallisto
    assert os.path.isdir(os.path.split(args.output)[0]), args.output
    min_cutoff = float(args.cutoff)

    # Parse expression
    expr = pd.read_table(args.kallisto)
    expr['target_id'] = expr['target_id'].str.replace(r'\.\d+', '', regex=True)
    expr = expr[['target_id', 'tpm']].set_index('target_id')
    id2expr = expr.to_dict()['tpm']
      
    # Parse neoantigen table
    var = pd.read_table(args.pvacseq)
    var['Transcript Expression'] = var['Transcript'].map(id2expr) # add tpm
    var = var[var['Transcript Expression'] >= min_cutoff] # filter by tpm
    var.to_csv(args.output, sep='\t', index=False) # save to output
