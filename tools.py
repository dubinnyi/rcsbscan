
import gzip
import sys, os
import warnings
import time
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from Bio.PDB import *
from Bio.Seq import Seq
from Bio.Data import IUPACData

import Bio.SVDSuperimposer

import multiprocessing as mp
LOCK = mp.Lock()

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
    # format argument does not work
    force_format = format = 'auto'
    if 'format' in kwargs.keys():
        force_format = kwargs['format']
    (file_id, filename_format)  = parse_structure_filename(filename)
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


def range_str_to_tuple(range_str):
    if not range_str or range_str == '*' or range_str == '.' or range_str.upper() == 'ALL':
        # Select all
        return None
    range_split = range_str.split('-')
    if len(range_split) == 0:
        # Select all
        return None
    elif len(range_split) == 1:
        res_start = int(range_split[0])
        res_finish = res_start
    else:
        res_start = int(range_split[0])
        res_finish = int(range_split[1])
    return (res_start, res_finish + 1)


def select_aa_residue_range_from_chain(chain, range_tuple):
    res_list_aa = chain_sequence_aaselect(chain)
    if range_tuple:
        res_selected = [res for res in res_list_aa
                        if res.get_id()[1] >= range_tuple[0] and res.get_id()[1] < range_tuple[1]]
        if not res_selected:
            raise IndexError("list index '{}' is out of range".format(range_tuple))
    else:
        res_selected = res_list_aa
    return res_selected


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
            if not (this_id == prev_id or this_id == prev_id + 1 ):
                return False
        prev_id = this_id
    return True

def get_res_and_atoms(chain_res_aa, resno, ref_res_list_len, ref_atom_names_set):
    try:
        res_list_to_fit = chain_res_aa[resno : resno + ref_res_list_len]
    except:
        return (None, None)
    if check_res_list_is_continuous(res_list_to_fit):
        atom_list_to_fit = select_atoms_from_res_list(res_list_to_fit, ref_atom_names_set)
        return (res_list_to_fit, atom_list_to_fit)
    else:
        return (None, None)


def recursive_expand(list_of_files):
    result = []
    for f in list_of_files:
        if os.path.isfile(f):
            result.append(f)
        elif os.path.isdir(f):
            subdirs = [os.path.join(f,d) for d in os.listdir(f)]
            result.extend(recursive_expand(subdirs))
    return result


