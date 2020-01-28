#!/usr/bin/python3

from tools import *


def main():
    arg_parser = argparse.ArgumentParser(
        usage='Print sequence information from pdb or mmcif file(s)')
    arg_parser.add_argument('struct', type=str, help='structures in pdb or mmcif format', nargs='+')
    arg_parser.add_argument('-p', '--print-header', action='store_true',
                            help='print header')
    arg_parser.add_argument('-c', '--print-chain', action='store_true',
                            help='print chain sequence in one-letter format')
    arg_parser.add_argument('-C', '--print-chain-full', action='store_true',
                            help='print full chain sequence (one line per residue)')
    arg_parser.add_argument('-w', '--pdb-warnings', action='store_true',
                            help='show structure parsing warnings', default=False)
    arg_parser.add_argument('-i', '--print-insertions', action='store_true',
                            help='print insertions in sequence')
    arg_parser.add_argument('--format', type=str, default='auto',
                            help='File format: pdb, cif or auto')

    args = arg_parser.parse_args()
    if not args.pdb_warnings:
        warnings.simplefilter("ignore", PDBConstructionWarning)

    if args.struct:
        for structf in args.struct:
            try:
                (file_id, format) = parse_structure_filename(structf)
                structure = get_structure_from_file(structf, format=args.format)
                if not structure:
                    continue
                if args.print_header:
                    for key, value in sorted(structure.header.items()):
                        print("{:20} : {}".format(key, value))
                for model in structure:
                    for chain in model:
                        # seq_ppb = chain_sequence_ppbuild(chain)
                        chain_aa = chain_sequence_aaselect(chain)
                        (seq_start, Seq_aa) = protein_one_letter_string(chain_aa)
                        res_ins = res_insertions(chain_aa)
                        seq_water = chain_water(chain)
                        # len_aa, len_ppb, len_wat, len_ins = 0, 0, 0, 0
                        # len_ppb = 0 if not seq_ppb else len(seq_ppb)
                        len_aa = 0 if not chain_aa else len(chain_aa)
                        len_aa_gaps = 0 if not chain_aa else len(Seq_aa)
                        len_ins = 0 if not res_ins else len(res_ins)
                        len_water = 0 if not seq_water else len(seq_water)
                        print("> {} model= {:2}, chain= {:1}, {:4} ({:4}) AA, start = {:4}, {:3} insertions, {:3} WAT".
                              format(file_id, model.get_id(), chain.get_id(),
                                     len_aa, len_aa_gaps, seq_start, len_ins, len_water))
                        if args.print_chain:
                            if chain_aa and len(chain_aa) > 0:
                                for i in range(0, len(Seq_aa), 60):
                                    print("   {}".format(Seq_aa[i:i + 60]))
                                    # if seq_ppb and len(seq_ppb)>0 :
                                    #    for i in range(0,len(seq_ppb),60):
                                    #        print("   {}".format(seq_ppb[i:i+60]))
                        if args.print_chain_full:
                            print("All residues in chain {}:".format(chain))
                            for residue in chain:
                                print(residue)
            except FileNotFoundError:
                eprint("File not found: {}".format(structf))


if __name__ == "__main__":
    main()
