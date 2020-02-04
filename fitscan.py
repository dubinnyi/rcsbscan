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
    sup = Superimposer()
    t0 = time.time()
    count_files = 0
    count_struct = 0
    count_res_tuple = 0
    if args.struct:
        for structf in args.struct:
            try:
                (file_id, format) = parse_structure_filename(structf)
                structure = get_structure_from_file(structf, format='auto')
                if not structure:
                    continue
                count_files = count_files + 1
                if args.print_header and args.verbose:
                    for key, value in sorted(structure.header.items()):
                        print("{:20} : {}".format(key, value))
                for model in structure:
                    for chain in model:
                        chain_res_aa = chain_sequence_aaselect(chain)
                        (seq_start, Seq_aa) = chain_one_letter_string(chain)
                        seq_water = chain_water(chain)
                        # len_aa, len_ppb, len_wat, len_ins = 0, 0, 0, 0
                        # len_ppb = 0 if not seq_ppb else len(seq_ppb)
                        len_aa = 0 if not chain_res_aa else len(chain_res_aa)
                        len_aa_gaps = 0 if not chain_res_aa else len(Seq_aa)
                        len_water = 0 if not seq_water else len(seq_water)
                        count_struct = count_struct + 1
                        if args.verbose:
                            print("> {} model= {:2}, chain= {:1}, {:4} ({:4}) AA, start = {:4}, {:3} WAT".
                                  format(file_id, model.get_id(), chain.get_id(),
                                         len_aa, len_aa_gaps, seq_start, len_water))
                            if len_aa:
                                for i in range(0, len(Seq_aa), 60):
                                    print("   {}".format(Seq_aa[i:i + 60]))
                        for resno in range(len_aa - s4fit.ref_res_list_len):
                            (res_to_fit, atoms_to_fit) = get_res_and_atoms(chain_res_aa, resno, s4fit.ref_res_list_len,
                                                                           s4fit.ref_atom_names_set)
                            if atoms_to_fit and len(atoms_to_fit) == len(s4fit.ref_atoms):
                                count_res_tuple = count_res_tuple + 1
                                sup.set_atoms(s4fit.ref_atoms, atoms_to_fit)
                                if sup.rms <= args.max_rms:
                                    (r_start, r_seq) = res_list_to_one_letter_string(res_to_fit)
                                    print("RMSD_HIT: {} model= {:>3}, chain= {:1} hit= {:>4}- {} -{:<4} rmsd= {:>6.4f}".
                                          format(file_id, model.get_id(), chain.get_id(),
                                                 r_start, r_seq, r_start + s4fit.ref_res_list_len - 1, sup.rms))

            except FileNotFoundError:
                eprint("Failed to read structure {}".format(structf))
                continue

    print("fitscan finished after {:>4d} sec.".format(int(time.time() - t0)))
    print("Statistics:")
    print("  {:>7} files".format(count_files))
    print("  {:>7} structures".format(count_struct))
    print("  {:>7} {}-length tuples of residues".format(count_res_tuple,s4fit.ref_res_list_len))


if __name__ == "__main__":
    main()
