import gzip
import sys, os
import warnings
import time
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import re

from Bio.PDB import *
from Bio.Seq import Seq
from Bio.Data import IUPACData

import Bio.SVDSuperimposer

import multiprocessing as mp

LOCK = mp.Lock()


class MPCounter(object):
    def __init__(self, initval=0):
        self.val = mp.Value('i', initval)
        self.lock = mp.Lock()

    def increment(self):
        with self.lock:
            self.val.value += 1

    def value(self):
        with self.lock:
            return self.val.value


def eprint(*args, **kwargs):
    with LOCK:
        print(*args, file=sys.stderr, **kwargs)


def mpprint(*args, **kwargs):
    with LOCK:
        print(*args, file=sys.stdout, **kwargs)


def opengz(filename, args: 'r'):
    if filename[-2:] == "gz":
        openfunc = gzip.open
        openarg = args + 't'
    else:
        openfunc = open
        openarg = args
    return openfunc(filename, openarg)


def parse_structure_filename(filename):
    id = os.path.basename(filename)
    format = 'auto'
    if id.startswith("pdb"):
        id = id[3:]
        format = 'pdb'
    if id.endswith(".gz"):
        id = id[0:-3]
    if id.endswith(".pdb"):
        id = id[0:-4]
        format = 'pdb'
    if id.endswith(".mmcif"):
        id = id[0:-6]
        format = 'cif'
    if id.endswith(".cif"):
        id = id[0:-4]
        format = 'cif'
    if id.endswith(".ent"):
        id = id[0:-4]
    id = id.upper()
    return (id, format)


def get_structure_from_file(filename, **kwargs):
    # format argument is unused
    force_format = format = 'auto'
    if 'format' in kwargs.keys():
        force_format = kwargs['format']
    (file_id, filename_format) = parse_structure_filename(filename)
    cif_parser = MMCIFParser()
    pdb_parser = PDBParser()
    if force_format != 'auto':
        format = force_format
    if format == 'auto':
        format = filename_format
    if format == 'pdb' or format == 'auto':
        parsers = [pdb_parser, cif_parser]
    else:
        parsers = [cif_parser, pdb_parser]
    for parser in parsers:
        try:
            with opengz(filename, 'r') as handle:
                structure = parser.get_structure(file_id, handle)
                if structure:
                    return structure
        except:
            continue
    eprint("File {} could not be read".format(filename))
    return None


def chain_sequence_ppbuild(chain):
    # Fails if chain is dicontinued
    ppb = PPBuilder()
    pp = ppb.build_peptides(chain)
    if pp and len(pp):
        return (pp[0].get_sequence())
    else:
        return (None)


def chain_sequence_aaselect(chain):
    # Ok if chain is discontinued
    res_list = chain.get_list()
    res_list_aa = [res for res in res_list if is_aa(res)]
    return res_list_aa


def res_oneletter(res):
    try:
        rname3 = res.get_resname()
        rname3 = rname3[0].upper() + rname3[1:3].lower()
        rname1 = IUPACData.protein_letters_3to1[rname3]
    except KeyError:
        rname1 = '.'
    return rname1


def res_list_to_one_letter_string(res_list):
    if not res_list:
        return (0, "")

    prev_id = None
    sequence_start = None
    seq = ''
    for r in res_list:
        this_id = r.get_id()[1]  # residue number, may repeat o have gaps
        if prev_id:
            if this_id == prev_id or this_id == prev_id + 1:
                seq += res_oneletter(r)
            elif this_id > prev_id + 1:
                for r1 in range(prev_id + 1, this_id):
                    seq += "."
                seq += res_oneletter(r)
        else:
            sequence_start = this_id
            seq += res_oneletter(r)
        prev_id = this_id
    return (sequence_start, Seq(seq))


def chain_one_letter_string(chain):
    res_list = chain_sequence_aaselect(chain)
    return res_list_to_one_letter_string(res_list)


def res_insertions(res_list):
    res_list_ins = [res for res in res_list if res.get_id()[2] != ' ']
    return res_list_ins


def chain_water(chain):
    # Ok if chain is discontinued
    res_list = chain.get_list()
    res_list_water = [res for res in res_list if res.get_id()[0] == 'W']
    return res_list_water


def chain_water_id(chain, water_id):
    # Ok if chain is discontinued
    res_list = chain.get_list()
    res_list_water = [res for res in res_list \
                      if res.get_id()[0] == 'W' and res.get_id()[1] == water_id]
    return res_list_water


def bio_str_to_resid(string, hetfield =' '):
    # Return residue in biopython-compatible format
    # insertion character after number is converted to icode
    # hetfield is ' ' for ordinary residues
    #
    # examples:
    # bio_str_to_resid('10')
    # (' ', 10, ' ')
    # bio_str_to_resid('11A')
    # (' ', 11, 'A')
    p = re.compile('(?P<resseq>\d+)(?P<icode>[A-Z]?)$')
    m = p.match(string)
    if not m:
        return None
    resseq = int(m.group('resseq'))
    icode = m.group('icode')
    icode = ' ' if not icode else icode
    return (hetfield, resseq, icode)


def bio_resid_to_str(resid):
    hetfield, resseq, icode = resid
    return "{}{}".format(resseq, icode.strip())

def range_str_to_tuple(range_str):
    if not range_str or range_str == '*' or range_str == '.' or range_str.upper() == 'ALL':
        # Select all
        return None
    range_split = range_str.split('-')
    if len(range_split) == 0:
        # Select all
        return None
    elif len(range_split) == 1:
        res_start = bio_str_to_resid(range_split[0])
        res_finish = res_start
    else:
        res_start = bio_str_to_resid(range_split[0])
        res_finish = bio_str_to_resid(range_split[1])
    #res_finish[1] += 1
    return (res_start, res_finish)


def select_residues_from_chain(chain, **kwargs):
    # chain: biopython chain
    # first_res: first residue id to select, e.g. (' ',1,' ')
    # last_res: last residue to select, e.g. (' ',10, 'A')
    # gap_count: select 'gap_count' residues before first_res or after last_res
    # select all residues if no argumants provided
    res_list_aa = chain_sequence_aaselect(chain)
    id_list = [res.get_id() for res in res_list_aa]

    gap_count = kwargs['gap_count'] if 'gap_count' in kwargs else 0
    first_res = kwargs['first_res'] if 'first_res' in kwargs else None
    last_res  = kwargs['last_res']  if 'last_res'  in kwargs else None

    first_index, last_index = 0, len(res_list_aa)
    if first_res:
        try:
            first_index = id_list.index(first_res)
        except:
            raise IndexError("Residue '{}' is out of sequence range".format(first_res))
    if last_res:
        try:
            last_index = id_list.index(last_res)
        except:
            raise IndexError("Residue '{}' is out of sequence range".format(last_res))

    if first_res and last_res:
        return res_list_aa[first_index: last_index + 1]
    elif first_res and gap_count:
        if first_index - gap_count < 0:
            gap_count = first_index
        return res_list_aa[first_index - gap_count : first_index]
    elif last_res and gap_count:
        return res_list_aa[last_index + 1 : last_index + 1 + gap_count]
    elif first_res and not gap_count:
        return res_list_aa[first_index :]
    elif last_res and not gap_count:
        return res_list_aa[: last_index + 1]
    else:
        return res_list_aa


def atom_str_to_set(atom_str):
    atom_set = set()
    if not atom_str or atom_str == '*' or atom_str == '.' or atom_str.upper() == 'ALL':
        return atom_set
    atom_list = atom_str.split(',')
    for a in atom_list:
        a = a.strip()
        atom_set.add(a)
    return atom_set


def select_atoms_from_res_list(res_list, ref_atom_names_set):
    atom_list = []
    for r in res_list:
        for a in r.get_list():
            if ref_atom_names_set and a.get_name() in ref_atom_names_set:
                atom_list.append(a)
            if not ref_atom_names_set:
                # set is empty, append all atoms
                atom_list.append(a)
    if ref_atom_names_set and not atom_list:
        raise IndexError("Atoms not present in structure")
    return atom_list


def check_res_list_is_continuous(res_list):
    prev_id = None
    for r in res_list:
        this_id = r.get_id()[1]  # residue number, may repeat o have gaps
        if prev_id:
            if not (this_id == prev_id or this_id == prev_id + 1):
                return False
        prev_id = this_id
    return True


def get_res_and_atoms(chain_res_aa, resno, ref_res_list_len, ref_atom_names_set):
    try:
        res_list_to_fit = chain_res_aa[resno: resno + ref_res_list_len]
    except:
        return (None, None)
    if check_res_list_is_continuous(res_list_to_fit):
        atom_list_to_fit = select_atoms_from_res_list(res_list_to_fit, ref_atom_names_set)
        return (res_list_to_fit, atom_list_to_fit)
    else:
        # eprint("Residue list is discontinous: {} - {}".format(res_list_to_fit[0].get_full_id(),
        #                                                      res_list_to_fit[-1].get_id()))
        return (None, None)


def recursive_expand(list_of_files_and_dirs):
    result = []
    for f in list_of_files_and_dirs:
        if os.path.isfile(f):
            result.append(f)
        elif os.path.isdir(f):
            subdirs = [os.path.join(f, d) for d in os.listdir(f)]
            result.extend(recursive_expand(subdirs))
    return result
