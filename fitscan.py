#!/usr/bin/python3 -u

from tools import *
from struct4fit import Struct4Fit
import argparse
import time

def main():
    arg_parser = argparse.ArgumentParser(
        usage='Fit reference structure to every substructute')
    arg_parser.add_argument('struct', type = str, help = 'structures in pdb or mmcif format', nargs = '+')
    arg_parser.add_argument('-p', '--print-header', action='store_true',
                            help='print header')
    arg_parser.add_argument('--ref-structure', type=str, help='Reference structure', required=True)
    arg_parser.add_argument('--ref-model', type=int, help='Model number in reference structure', default = None)
    arg_parser.add_argument('--ref-chain', type=str, help='Chain in reference structure', default = None)
    arg_parser.add_argument('--ref-residues', type=str, help='Residue range in reference structure', default = None)
    arg_parser.add_argument('--ref-atoms', type=str, help='Atoms in reference structure',
                            default='N,CA,C,O')
    arg_parser.add_argument('-w', '--pdb-warnings', action='store_true',
                            help='show structure parsing warnings', default=False)
    arg_parser.add_argument('-v', '--verbose', action='store_true',
                            help='Verbose output', default=False)
    arg_parser.add_argument('--max-rms', type=float, help='Maximum RMSD to print, ingore otherwise', default=1.0)
    args = arg_parser.parse_args()

    if not args.pdb_warnings:
        warnings.simplefilter("ignore", PDBConstructionWarning)

    #ref_struct=get_structure_from_file(args.ref_structure)

    s4fit = Struct4Fit(args.ref_structure, args.ref_model, args.ref_chain,
                       args.ref_residues,  args.ref_atoms, args.verbose)
    if not s4fit:
        quit(-1)

    ##############################
    #  Start fit scan
    ##############################

    print("Start fit scan")

    s4fit.reset_timer()
    if args.struct:
        for structf in args.struct:
            try:
                s4fit.scan(structf, args.max_rms)
            except FileNotFoundError:
                eprint("Failed to read structure {}".format(structf))
                continue

    s4fit.print_stats()


if __name__ == "__main__":
    main()
