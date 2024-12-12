#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse
import sys
from pathlib import Path
import glob




def parse_args(argv=None):
    """Define and immediately parse command line arguments."""

    parser = argparse.ArgumentParser(
        description="create a samplesheet from smart-seq samplesheet of EI format",
    )
    parser.add_argument(
        "--input",
        metavar="FILE_SAMPLESHEET",
        type=Path,
        help="Input smart-seq2 samplesheet file.",
        required=True,
    )
    parser.add_argument(
        "--seqdir",
        default='.',
        help="Directory of raw smart-seq2 data.", # contains subfolers of each cell sample
    )
    parser.add_argument(
        "--parity",
        default='paired',
        choices=['paired', 'single'],
        help="Choose sinlge-end reads or paired-end reads",
    )
    parser.add_argument(
        "--group",
        default='meta_1',
        help="Choose a meta column for group in final samplesheet", # contains subfolers of each cell sample
    )    
    parser.add_argument(
        "--outfile",
        metavar="OUTFILE",
        type=Path,
        help="Output batch file.",
        required=True,
    )                    
    return parser.parse_args(argv)


def main(argv=None):

    args = parse_args(argv)

    if not args.input.is_file():
        print(f"The given input file {args.input} was not found!")
        sys.exit(2)

    samplesheet = pd.read_csv(args.input)
    samplesheet.columns = samplesheet.columns.str.lower()
    samplesheet = samplesheet[samplesheet['control'] == False]    

    # samples_to_remove = []

    with open(args.outfile, 'w') as fh:
        if args.parity == 'single':
            fh.write(','.join(['sample', 'fastq_1', 'plate', 'group']) + '\n')
        elif args.parity == 'paired':
            fh.write(','.join(['sample', 'fastq_1', 'fastq_2', 'plate', 'group']) + '\n')
        for i, sample in enumerate(samplesheet['sample_id']):
            plate = samplesheet['sample_plate'].iloc[i]
            group = samplesheet[args.group].iloc[i]            
            if args.parity == 'single':
                reads = glob.glob(str(Path(args.seqdir, sample,  '*.fastq.gz')))
                if reads:
                    fh.write(','.join([sample, reads[0], plate, group]) + '\n')
                # else:
                #     samples_to_remove += [sample]
            elif args.parity == 'paired':
                read1 = glob.glob(str(Path(args.seqdir, sample,  '*_R1.fastq.gz')))
                read2 = glob.glob(str(Path(args.seqdir, sample,  '*_R2.fastq.gz')))
                if read1 and read2:
                    fh.write(','.join([sample, read1[0], read2[0], plate, group]) + '\n')
                # else:
                #     samples_to_remove += [sample]
    
    # samplesheet = samplesheet[~samplesheet["sample_id"].isin(samples_to_remove)]




if __name__ == "__main__":
    sys.exit(main())
