#!/usr/bin/python3

from tools import *
import argparse

def main():
    arg_parser = argparse.ArgumentParser(
        usage='Fit reference structure to every substructute')
    arg_parser.add_argument('struct', type = str, help = 'structures in pdb or mmcif format', nargs = '+')
    arg_parser.add_argument('-p', '--print-header', action='store_true',
                            help='print header')
    arg_parser.add_argument('--ref-structure', type=str, help='Reference structure', required=True)
    arg_parser.add_argument('--ref-chain', type=str, help='Chain in reference structure')
    arg_parser.add_argument('--ref-model', type=int, help='Model number in reference structure')
    arg_parser.add_argument('--ref-residues', type=str, help='Residue range in reference structure')
    arg_parser.add_argument('--ref-atoms', type=str, help='Atoms in reference structure',
                            default='N,CA,C,O')
    arg_parser.add_argument('-w', '--pdb-warnings', action='store_true',
                            help='show structure parsing warnings', default=False)
    arg_parser.add_argument('--max-rms', type=float, help='Maximum RMSD to print, ingore otherwise', default=1.0)
    args = arg_parser.parse_args()

    if not args.pdb_warnings:
        warnings.simplefilter("ignore", PDBConstructionWarning)

    ref_struct=get_structure_from_file(args.ref_structure)

    # Select model
    if args.ref_model:
        ref_mdl = ref_struct[args.ref_model]
    else:
        ref_mdl = ref_struct[0]

    # Select chain
    if args.ref_chain:
        ref_ch_name = args.ref_chain
    else:
        chain_id_list = [ch.get_id() for ch in ref_mdl.get_list()]
        ref_ch_name = chain_id_list[0]
        print("Automatically select first chain: \'{}\'".format(ref_ch_name))
    ref_ch = ref_mdl[ref_ch_name]

    (seq_start, Seq_aa) = chain_one_letter_string(ref_ch)
    print("{} residues in chain \'{}\', first residue is {}".format(len(Seq_aa), ref_ch_name, seq_start))
    for i in range(0, len(Seq_aa), 60):
        print("   {}".format(Seq_aa[i:i + 60]))

    if args.ref_residues:
        res_range_tuple = range_str_to_tuple(args.ref_residues)
    else:
        res_range_tuple = None # Select all
    ref_res_list = select_aa_residue_range_from_chain(ref_ch, res_range_tuple)
    ref_res_list_len= len(ref_res_list)

    (seq2_start, Seq2_aa) = res_list_to_one_letter_string(ref_res_list)
    print("{} residues selected from chain \'{}\'".format(len(Seq2_aa), ref_ch_name))
    for i in range(0, len(Seq2_aa), 60):
        print("   {}".format(Seq2_aa[i:i + 60]))

    if args.ref_atoms:
        ref_select_atoms_str = args.ref_atoms
        ref_atom_names_set = atom_str_to_set(ref_select_atoms_str)
    else:
        ref_select_atoms_str = ''
        ref_atom_names_set = set()
    ref_atoms = select_atoms_from_res_list(ref_res_list, ref_atom_names_set)
    print("{} atoms selected by filter \'{}\'".format(len(ref_atoms), ref_select_atoms_str))
    for a in ref_atoms:
        print(a.get_full_id())

    if not ref_atoms:
        eprint("No atoms selected. Exiting")
        return(-1)

    ##############################
    #  Start fit scan
    ##############################

    print("Start fit scan")
    sup = Superimposer()
    if args.struct:
        for structf in args.struct:
            try:
                (file_id, format) = parse_structure_filename(structf)
                structure = get_structure_from_file(structf, format='auto')
                if not structure:
                    continue
                if args.print_header:
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
                        print("> {} model= {:2}, chain= {:1}, {:4} ({:4}) AA, start = {:4}, {:3} WAT".
                              format(file_id, model.get_id(), chain.get_id(),
                                     len_aa, len_aa_gaps, seq_start, len_water))
                        if len_aa:
                            for i in range(0, len(Seq_aa), 60):
                                print("   {}".format(Seq_aa[i:i + 60]))
                        for resno in range(len_aa - ref_res_list_len):
                            (res_to_fit, atoms_to_fit) = get_res_and_atoms(chain_res_aa, resno, ref_res_list_len,
                                                             ref_atom_names_set)
                            if atoms_to_fit and len(atoms_to_fit) == len(ref_atoms):
                                sup.set_atoms(ref_atoms, atoms_to_fit)
                                if sup.rms <= args.max_rms:
                                    (r_start, r_seq) = res_list_to_one_letter_string(res_to_fit)
                                    print("RMSD HIT: {} model={:2}, chain={:1} {:4>}-{} rmsd = {:>6.4f}".
                                          format(file_id, model.get_id(), chain.get_id(), r_start, r_seq, sup.rms))

            except FileNotFoundError:
                eprint("Failed to read structure {}".format(structf))
                continue


if __name__ == "__main__":
    main()
