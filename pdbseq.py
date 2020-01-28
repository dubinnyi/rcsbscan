#!/usr/bin/python3

import argparse
import gzip
import sys, os
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from Bio.PDB import *
from Bio.Seq import Seq
from Bio.Data import IUPACData


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


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
    structure = None
    for parser in parsers:
        try:
            with opengz(filename, 'r') as handle:
                structure = parser.get_structure(file_id, handle)
                if structure:
                    return structure
        except:
            continue
    return structure


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


def protein_one_letter_string(res_list):
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


def res_insertions(res_list):
    res_list_ins = [res for res in res_list if res.get_id()[2] != ' ']
    return res_list_ins


def chain_water(chain):
    # Ok if chain is discontinued
    res_list = chain.get_list()
    res_list_water = [res for res in res_list if res.get_id()[0] == 'W']
    return res_list_water


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
