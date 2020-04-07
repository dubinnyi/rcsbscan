#!/usr/bin/python3

from scan.biotools import *
import argparse

from Bio.PDB.StructureBuilder import StructureBuilder


def pdb_extract(structure, extract_model, extract_chain, res_range_tuple, water_id):

    structure_builder = StructureBuilder()
    structure_builder.init_structure('pdb_extract')
    structure_builder.set_line_counter(0)
    start_resseq = 0
    current_model_id = 0
    line_counter = 0

    for model in structure:
        if extract_model and model.get_id() != extract_model:
            continue

        structure_builder.init_model(current_model_id, current_model_id)
        current_model_id += 1

        for chain in model:
            if extract_chain and chain.get_id() != extract_chain:
                continue
            structure_builder.init_seg(' ')
            structure_builder.init_chain(chain.get_id())
            new_resseq = start_resseq
            for residue in chain:
                hetero_flag, resseq, icode = residue.get_id()
                if res_range_tuple and ( not hetero_flag or hetero_flag == ' ' ) :
                    if res_range_tuple[0] > resseq or res_range_tuple[1] <= resseq :
                        continue
                if hetero_flag == 'W':
                    if not water_id:
                        continue
                    elif resseq != water_id:
                        continue
                if hetero_flag.startswith('H'):
                    continue
                new_resseq += 1
                structure_builder.init_residue(residue.get_resname(), hetero_flag, new_resseq, icode)
                for atom in residue:
                    structure_builder.init_atom(atom.get_name(), atom.get_coord(),
                                                atom.get_bfactor(), atom.get_occupancy(), atom.get_altloc(),
                                                atom.get_fullname())
                    structure_builder.set_line_counter(line_counter)
                    line_counter += 1
    out_structure = structure_builder.get_structure()
    return out_structure


def main():
    arg_parser = argparse.ArgumentParser(
        usage='Print sequence information from pdb or mmcif file(s)')
    arg_parser.add_argument('struct', type=str, help='Structure in pdb or mmcif format')
    arg_parser.add_argument('out', type=str, help='Output file for extructed substructure')
    arg_parser.add_argument('--model', type=int, help='Model number to extract from structure', default=None)
    arg_parser.add_argument('--chain', type=str, help='Chain in structure', default=None)
    arg_parser.add_argument('--residues', type=str, help='Residue range in structure', default = None)
    arg_parser.add_argument('--water', type=int, help='Water molecule to extract', default=None)
    arg_parser.add_argument('-w', '--pdb-warnings', action='store_true',
                            help='show structure parsing warnings', default=False)

    args = arg_parser.parse_args()

    if not args.pdb_warnings:
        warnings.simplefilter("ignore", PDBConstructionWarning)

    res_range_tuple = None
    if args.residues:
        res_range_tuple = range_str_to_tuple(args.residues)

    try:
        structure = get_structure_from_file(args.struct, format='pdb')
        out_structure = pdb_extract(structure, args.model, args.chain, res_range_tuple, args.water)
    except FileNotFoundError:
        eprint("File not found: {}".format(args.struct))
        return None

    try:
        out_pdb = args.out
        pdbio = PDBIO(True)
        pdbio.set_structure(out_structure)
        pdbio.save(out_pdb, write_end=True)
    except:
        eprint("Could not write to file {}".format(out_pdb))


if __name__ == "__main__":
    main()
