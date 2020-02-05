#!/usr/bin/python3 -u

from tools import *

class FitCounter:
    def __init__(self, name):
        self.name = name
        self.counters_reset()
        self.timer_reset()

    def timer_reset(self):
        self.t0 = None
        self.t1 = None
        self.total_time = 0

    def timer_start(self):
        self.t0 = time.time()
        self.t1 = None
        self.total_time = 0

    def timer_finish(self):
        self.t1 = time.time()
        self.total_time = self.total_time + (self.t1 - self.t0)

    def counters_reset(self):
        self.count_files = 0
        self.count_struct = 0
        self.count_res_tuple = 0
        self.count_hits = 0
        self.timer_start()

    def new_file(self):
        self.count_files = self.count_files + 1

    def new_struct(self):
        self.count_struct = self.count_struct + 1

    def new_res_tuple(self):
        self.count_res_tuple = self.count_res_tuple + 1

    def new_hit(self):
        self.count_hits = self.count_hits + 1

    def __str__(self):
        out = []
        out.append("fitscan statistics after {:>4d} sec:".format(int(self.total_time)))
        out.append("  {:>7} files".format(self.count_files))
        out.append("  {:>7} structures".format(self.count_struct))
        out.append("  {:>7} hits".format(self.count_hits))
        out.append("  {:>7} {}".format(self.count_res_tuple, self.name))
        return "\n".join(out)

    def print_stats(self):
        self.timer_finish()
        print(str(self))

    def __add__(self, counter):
        self.total_time =      self.total_time +      counter.total_time
        self.count_files =     self.count_files +     counter.count_files
        self.count_struct =    self.count_struct +    counter.count_struct
        self.count_res_tuple = self.count_res_tuple + counter.count_res_tuple
        self.count_hits =      self.count_hits +      counter.count_hits
        return self



class Struct4Fit:
    def __init__(self, struct, model, chain, residues, atoms, verbose):
        self.ok_flag = False
        self.verbose = verbose
        self.struct = get_structure_from_file(struct)
        self.sup = Superimposer()
        if self.struct:

            # Select model
            if model:
                self.model = model
                self.ref_mdl = self.struct[model]
            else:
                self.model = 0
                self.ref_mdl = self.struct[0]

            # Select chain
            if chain:
                self.ref_ch_name = chain
            else:
                chain_id_list = [ch.get_id() for ch in self.ref_mdl.get_list()]
                self.ref_ch_name = chain_id_list[0]
                if self.verbose:
                    print("Automatically select first chain: \'{}\'".format(self.ref_ch_name))
            self.ref_ch = self.ref_mdl[self.ref_ch_name]

            (seq_start, Seq_aa) = chain_one_letter_string(self.ref_ch)
            if self.verbose:
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
            self.ref_seq = Seq2_aa
            self.ref_seq_start = seq2_start
            self.ref_seq_len = len(self.ref_seq)
            self.ref_seq_end = self.ref_seq_start + self.ref_seq_len - 1
            self.ref_seq_formatted_60 = []
            for i in range(0, self.ref_seq_len, 60):
                self.ref_seq_formatted_60.append("{}".format(self.ref_seq[i:i + 60]))
            if self.verbose:
                print("{} residues selected from chain \'{}\'".format(self.ref_seq_len, self.ref_ch_name))
                print("\n".join(self.ref_seq_formatted_60))

            if atoms:
                self.ref_select_atoms_str = atoms
                self.ref_atom_names_set = atom_str_to_set(self.ref_select_atoms_str)
            else:
                self.ref_select_atoms_str = ''
                self.ref_atom_names_set = set()
            self.ref_atoms = select_atoms_from_res_list(self.ref_res_list, self.ref_atom_names_set)
            if self.verbose:
                print("{} atoms selected by filter \'{}\'".format(len(self.ref_atoms), self.ref_select_atoms_str))
                for a in self.ref_atoms:
                    print(a.get_full_id())

            if not self.ref_atoms:
                eprint("No atoms selected in reference structure. Exiting")
            else:
                self.ok_flag = True
                print("REF_4FIT: {} model= {:>3}, chain= {:1} seq= {:>4} {} {:<4} atoms={} fit_atoms={}".
                      format(self.struct.get_id(), self.model, self.ref_ch_name,
                      self.ref_seq_start, self.ref_seq, self.ref_seq_end, self.ref_select_atoms_str,
                             len(self.ref_atoms)))
        self.counter = FitCounter("{}-length tuples of residues".format(self.ref_res_list_len))
        self.counter.counters_reset()


    def scan(self, structf, max_rms):
        (file_id, format) = parse_structure_filename(structf)
        structure = get_structure_from_file(structf, format='auto')
        if not structure:
            return
        self.counter.new_file()
        if self.verbose:
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
                self.counter.new_struct()
                if self.verbose:
                    print("> {} model= {:2}, chain= {:1}, {:4} ({:4}) AA, start = {:4}, {:3} WAT".
                          format(file_id, model.get_id(), chain.get_id(),
                                 len_aa, len_aa_gaps, seq_start, len_water))
                    if len_aa:
                        for i in range(0, len(Seq_aa), 60):
                            print("   {}".format(Seq_aa[i:i + 60]))
                for resno in range(len_aa - self.ref_res_list_len):
                    (res_to_fit, atoms_to_fit) = get_res_and_atoms(chain_res_aa, resno, self.ref_res_list_len,
                                                                   self.ref_atom_names_set)
                    if atoms_to_fit and len(atoms_to_fit) == len(self.ref_atoms):
                        self.counter.new_res_tuple()
                        self.sup.set_atoms(self.ref_atoms, atoms_to_fit)
                        if self.sup.rms <= max_rms:
                            self.counter.new_hit()
                            (r_start, r_seq) = res_list_to_one_letter_string(res_to_fit)
                            print("RMSD_HIT: {} model= {:>3}, chain= {:1} hit= {:>4} {} {:<4} rmsd= {:>6.4f}".
                                  format(file_id, model.get_id(), chain.get_id(),
                                         r_start, r_seq, r_start + self.ref_res_list_len - 1, self.sup.rms))

    def __bool__(self):
        return(self.ok_flag)
