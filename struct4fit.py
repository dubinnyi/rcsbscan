#!/usr/bin/python3 -u

from tools import *

class Struct4Fit:
    def __init__(self, struct, model, chain, residues, atoms, verbose):
        self.ok_flag = False
        self.struct = get_structure_from_file(struct)
        if self.struct:

            # Select model
            if model:
                self.ref_mdl = self.struct[model]
            else:
                self.ref_mdl = self.struct[0]

            # Select chain
            if chain:
                self.ref_ch_name = chain
            else:
                chain_id_list = [ch.get_id() for ch in self.ref_mdl.get_list()]
                self.ref_ch_name = chain_id_list[0]
                print("Automatically select first chain: \'{}\'".format(self.ref_ch_name))
            self.ref_ch = self.ref_mdl[self.ref_ch_name]

            (seq_start, Seq_aa) = chain_one_letter_string(self.ref_ch)
            print("{} residues in chain \'{}\', first residue is {}".format(len(Seq_aa), self.ref_ch_name, seq_start))
            for i in range(0, len(Seq_aa), 60):
                print("   {}".format(Seq_aa[i:i + 60]))

            if residues:
                res_range_tuple = range_str_to_tuple(residues)
            else:
                res_range_tuple = None  # Select all
            self.ref_res_list = select_aa_residue_range_from_chain(self.ref_ch, res_range_tuple)
            self.ref_res_list_len = len(self.ref_res_list)

            (seq2_start, Seq2_aa) = res_list_to_one_letter_string(self.ref_res_list)
            print("{} residues selected from chain \'{}\'".format(len(Seq2_aa), self.ref_ch_name))
            for i in range(0, len(Seq2_aa), 60):
                print("   {}".format(Seq2_aa[i:i + 60]))

            if atoms:
                self.ref_select_atoms_str = atoms
                self.ref_atom_names_set = atom_str_to_set(self.ref_select_atoms_str)
            else:
                ref_select_atoms_str = ''
                self.ref_atom_names_set = set()
            self.ref_atoms = select_atoms_from_res_list(self.ref_res_list, self.ref_atom_names_set)
            print("{} atoms selected by filter \'{}\'".format(len(self.ref_atoms), self.ref_select_atoms_str))
            if verbose:
                for a in self.ref_atoms:
                    print(a.get_full_id())

            if not self.ref_atoms:
                eprint("No atoms selected in reference structure. Exiting")
            else:
                self.ok_flag = True

    def __bool__(self):
        return(self.ok_flag)
