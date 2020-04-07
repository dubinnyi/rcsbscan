#!/usr/bin/python3 -u

from scan.biotools import *
from collections import namedtuple
from pdbextract import pdb_extract


class RMSD_Hit(namedtuple('RMSD_Hit', ['pdb', 'chain', 'model', 'res_start', 'res_end', 'hit_sequence', 'rmsd'])):
    def __str__(self):
        out = "RMSD_HIT: {} model= {:>3}, chain= {:1} hit= {:>4} {} {:<4} rms= {:>6.4f}"
        out = out.format(self.pdb, self.model, self.chain,
                         self.res_start, self.hit_sequence, self.res_end, self.rmsd)
        return out


class Water_Hit(namedtuple('Water_Hit', ['water_id', 'water_rms'])):
    def __init__(self):
        self.water_id = None
        self.water_rms = None

    def __str__(self):
        str = ""
        if self.water_id and self.water_rms:
            str = "WAT= {:4} wat_rms= {:>6.4f}"
            str.format(format(self.water_id, self.water_rms))
        return str


class Struct4Fit:

    def __init__(self, struct_file, model, chain, residues, atoms, verbose, max_rms):
        self.ok_flag = False
        self.verbose = verbose
        self.struct = get_structure_from_file(struct_file)
        self.sup = Superimposer()
        self.max_rms = max_rms

        self.water_id = None
        self.water_res = []
        self.water_atoms = []
        self.water_vector = None

        self.out_filename = None
        self.out_filehandle = None
        self.pdbio = None
        self.Xray_max_resolution = None
        self.Xray_only = None

        if self.struct:

            # Select model
            try:
                if model:
                    self.model = model
                else:
                    self.model = 0
                self.ref_mdl = self.struct[self.model]
            except:
                eprint("Could not get model '{}' from file '{}'".format(model, struct_file))
                return

            # Select chain
            try:
                if chain:
                    self.ref_ch_name = chain
                else:
                    chain_id_list = [ch.get_id() for ch in self.ref_mdl.get_list()]
                    self.ref_ch_name = chain_id_list[0]
                    if self.verbose:
                        print("Automatically select first chain: \'{}\'".format(self.ref_ch_name))
                self.ref_ch = self.ref_mdl[self.ref_ch_name]
            except:
                eprint("Could not get chain '{}' from model '{}' in file '{}'".format(chain, self.model, struct_file))
                return

            (seq_start, Seq_aa) = chain_one_letter_string(self.ref_ch)
            if self.verbose:
                print(
                    "{} residues in chain \'{}\', first residue is {}".format(len(Seq_aa), self.ref_ch_name, seq_start))
                for i in range(0, len(Seq_aa), 60):
                    print("   {}".format(Seq_aa[i:i + 60]))

            # Select residues
            try:
                if residues:
                    self.res_range_tuple = range_str_to_tuple(residues)
                else:
                    self.res_range_tuple = None  # Select all
                if self.res_range_tuple:
                    self.res_range_str = "{}-{}".format(self.res_range_tuple[0], self.res_range_tuple[1] - 1)
                else:
                    self.res_range_str = "All"
                self.ref_res_list = select_aa_residue_range_from_chain(self.ref_ch, self.res_range_tuple)
                self.ref_res_list_len = len(self.ref_res_list)
            except:
                eprint("Could not get residues '{}' from chain '{}' of model '{}' in file '{}'".
                       format(residues, self.ref_ch_name, self.model, struct_file))
                return

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

            # Select atoms
            try:
                if atoms:
                    self.ref_select_atoms_str = atoms
                    self.ref_atom_names_set = atom_str_to_set(self.ref_select_atoms_str)
                else:
                    self.ref_select_atoms_str = '*'
                    self.ref_atom_names_set = set()  # Select all atoms
                self.ref_atoms = select_atoms_from_res_list(self.ref_res_list, self.ref_atom_names_set)
                if self.verbose:
                    print("{} atoms selected by filter \'{}\'".format(len(self.ref_atoms), self.ref_select_atoms_str))
                    for a in self.ref_atoms:
                        print(a.get_full_id())
            except:
                eprint("No atoms selected in reference structure:")
                eprint("Could not find atoms '{}' in residues '{}' from chain '{}' of model '{}' in file '{}'".
                       format(self.ref_select_atoms_str, self.res_range_str, self.ref_ch_name, self.model, struct_file))
                return

            # Success: some atoms selected
            self.ok_flag = True
            self.info = "REF_4FIT: {} model= {:>3}, chain= {:1} seq= {:>4} {} {:<4} "
            self.info = self.info + "atoms={} fit_atoms={} max_rms={:>6.4f}"
            self.info = self.info.format(self.struct.get_id(), self.model, self.ref_ch_name,
                                         self.ref_seq_start, self.ref_seq, self.ref_seq_end, self.ref_select_atoms_str,
                                         len(self.ref_atoms), self.max_rms)

            # print(self.info)
            self.result_name = "tuples of {} residues superimposed and rms of atoms {} evaluated".format(
                self.ref_res_list_len, self.ref_select_atoms_str, self.max_rms)

    def set_water(self, water_id, water_max_rms):
        self.water_id = None if not water_id else water_id
        self.water_max_rms = water_max_rms

        # Select water
        try:
            if self.water_id:
                # TODO: Support for multiple water molecules (?)
                water_reslist = chain_water(self.ref_ch)
                res_water_id = int(self.water_id)
                water_res = [res for res in water_reslist if res.get_id()[1] == res_water_id]
                if water_res:
                    self.water_res = [water_res[0]]
                    self.water_atoms_set = set('O')
                    self.water_atoms = select_atoms_from_res_list(self.water_res, self.water_atoms_set)
                    if self.water_atoms:
                        self.water_vector = self.water_atoms[0].get_vector()
                        print("Water {} selected by {}".format(self.water_res[0].get_id(), self.water_id))
                        print("Water coordinates: {}".format(self.water_vector))
                    else:
                        eprint("Coorinated of water molecule '{}' not found".format(self.water_id))
                else:
                    eprint("Water molecule '{}' not found".format(self.water_id))
        except:
            eprint("Can't select water molecule {}".format(self.water_id))
            return

        if self.water_vector:
            self.info += " WAT= {:<4}".format(self.water_id)

    def set_write(self, filename):
        self.out_filename = filename
        out_filehandle = open(self.out_filename, 'w')
        out_filehandle.write('')
        out_filehandle.close()
        self.info += "\nSaving PDB HITS to \'{}\'".format(self.out_filename)
        # self.pdbio = PDBIO(True)

    def set_Xray(self, max_resolution, xray_only):
        self.Xray_max_resolution = max_resolution
        self.Xray_only = xray_only
        self.info += "\nScan X-ray structures with resolution less then {}". \
            format(self.Xray_max_resolution)
        if self.Xray_only:
            self.info += "\nOnly X-RAY structures will be scanned"
        else:
            self.info += "\nAll non X-RAY structures will be scanned regrdless of resolution"

    def __str__(self):
        return self.info

    def check_method_and_resolution(self, structure):
        file_id = structure.get_id()
        xray_resolution = None
        if self.Xray_max_resolution and structure.header and 'resolution' in structure.header:
            xray_resolution = structure.header['resolution']
            if xray_resolution and self.Xray_max_resolution < xray_resolution:
                mpprint("Skip      {}: X-ray resolution {:4.2f} > {:4.2f}".
                        format(file_id, xray_resolution, self.Xray_max_resolution))
                return False
            elif not xray_resolution and not self.Xray_only:
                if self.verbose:
                    mpprint("Scan      {}: pdb does not have resolution. Processing.".
                            format(file_id))
                    return True
        if self.Xray_only and not xray_resolution:
            mpprint("Skip      {}: No resolution, Xray_only = True".
                    format(file_id))
            return False
        return True

    def scan(self, structf, counter):
        counter.name = "{}-length tuples of residues".format(self.ref_res_list_len)
        counter.new_file()
        try:
            structure = get_structure_from_file(structf, format='auto')
            file_id = structure.get_id()
            if not structure:
                counter.new_error()
                return

            if not self.check_method_and_resolution(structure):
                counter.new_skipped()
                return

            if self.verbose:
                for key, value in sorted(structure.header.items()):
                    print("{:20} : {}".format(key, value))
            for model in structure:
                for chain in model:
                    chain_res_aa = chain_sequence_aaselect(chain)
                    (seq_start, Seq_aa) = chain_one_letter_string(chain)
                    seq_water = chain_water(chain)
                    len_aa = 0 if not chain_res_aa else len(chain_res_aa)
                    len_aa_gaps = 0 if not chain_res_aa else len(Seq_aa)
                    len_water = 0 if not seq_water else len(seq_water)
                    counter.new_struct()
                    if self.verbose:
                        mpprint("> {} model= {:2}, chain= {:1}, {:4} ({:4}) AA, start = {:4}, {:3} WAT".
                                format(file_id, model.get_id(), chain.get_id(),
                                       len_aa, len_aa_gaps, seq_start, len_water))
                        if len_aa:
                            for i in range(0, len(Seq_aa), 60):
                                mpprint("   {}".format(Seq_aa[i:i + 60]))
                    if len_aa < self.ref_res_list_len:
                        continue
                    for resno in range(len_aa - self.ref_res_list_len):
                        (res_to_fit, atoms_to_fit) = get_res_and_atoms(chain_res_aa, resno, self.ref_res_list_len,
                                                                       self.ref_atom_names_set)
                        if atoms_to_fit and len(atoms_to_fit) == len(self.ref_atoms):
                            counter.new_res_tuple()
                            # print(atoms_to_fit)
                            self.sup.set_atoms(self.ref_atoms, atoms_to_fit)
                            if self.sup.rms <= self.max_rms:
                                # apply rotation to all atoms in pdb hit
                                # All atoms in residues, not only N,CA,C,O
                                all_hit_atoms = select_atoms_from_res_list(res_to_fit, set())
                                ###################################
                                self.sup.apply(all_hit_atoms)
                                ##################################
                                # check water match
                                water_match_str = ""
                                water_match_id = None
                                if self.water_vector:
                                    water_atoms = select_atoms_from_res_list(seq_water, self.water_atoms_set)
                                    if water_atoms:
                                        # water_ids = [a.full_id() for a in water_atoms]
                                        self.sup.apply(water_atoms)
                                        # print("Water molecules: {}".format(len(water_atoms)))
                                        # print("Water molecules tranformed: {}".format(len(water_atoms)))
                                        water_match = [water_atom for water_atom in water_atoms \
                                                       if water_atom - self.water_atoms[0] < self.water_max_rms]
                                        if water_match:
                                            water_match_0 = water_match[0]
                                            water_match_residue = water_match_0.get_parent()
                                            water_match_id = water_match_residue.get_id()[1]
                                            # water_match_id = water_match_id[1]
                                            water_rms = water_match[0] - self.water_atoms[0]
                                            water_match_str = "WAT= {:4} wat_rms= {:>6.4f}".format(water_match_id,
                                                                                                   water_rms)
                                counter.new_hit()
                                (r_start, r_seq) = res_list_to_one_letter_string(res_to_fit)
                                hit = RMSD_Hit(file_id, chain.get_id(), model.get_id(),
                                               r_start, r_start + self.ref_res_list_len - 1, str(r_seq), self.sup.rms)
                                mpprint("{} {}".format(hit, water_match_str))
                                if self.out_filename:
                                    with LOCK:
                                        with open(self.out_filename, 'a') as out_pdb:
                                            hit_res_tuple_start = res_to_fit[0].get_id()[1]
                                            hit_res_tuple = (
                                                hit_res_tuple_start, hit_res_tuple_start + self.ref_res_list_len)
                                            out_structure = pdb_extract(structure, model.get_id(), chain.get_id(),
                                                                        hit_res_tuple, water_match_id)
                                            out_pdb.write("REMARK 777 {} {}\n".format(hit, water_match_str))
                                            pdbio = PDBIO(True)
                                            pdbio.set_structure(out_structure)
                                            pdbio.save(out_pdb, write_end=False)
                                            out_pdb.flush()
                                            out_pdb.close()
        except IndexError as e:
            if self.verbose:
                eprint("Error in struct= {}, model= {}, chain= {}: {}".format(
                    file_id, model.get_id(), chain.get_id(), e))
                counter.new_error()

    def __bool__(self):
        return (self.ok_flag)
