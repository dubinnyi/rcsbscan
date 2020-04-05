#!/usr/bin/python3

from tools import *
import argparse

from Bio.PDB.StructureBuilder import StructureBuilder

def main():
    arg_parser = argparse.ArgumentParser(
        usage='Print sequence information from pdb or mmcif file(s)')
    arg_parser.add_argument('struct', type=str, help='Structure in pdb or mmcif format')
    arg_parser.add_argument('out', type=str, help='Output file for extructed substructure')
    arg_parser.add_argument('--model', type=int, help='Model number to extract from structure', default=None)
    arg_parser.add_argument('--chain', type=str, help='Chain in structure', default=None)
    arg_parser.add_argument('--residues', type=str, help='Residue range in structure', default = None)
    arg_parser.add_argument('-w', '--pdb-warnings', action='store_true',
                            help='show structure parsing warnings', default=False)

    args = arg_parser.parse_args()

    if not args.pdb_warnings:
        warnings.simplefilter("ignore", PDBConstructionWarning)

    res_range_tuple = None
    if args.residues:
        res_range_tuple = range_str_to_tuple(args.residues)

    structf = args.struct
    structure_builder = StructureBuilder()
    structure_builder.init_structure(args.out)
    structure_builder.set_line_counter(0)
    start_resseq = 0
    try:
        (file_id, format) = parse_structure_filename(structf)
        structure = get_structure_from_file(structf, format='pdb')

        current_model_id = 0
        line_counter = 0
        for model in structure:
            if  args.model and model != args.model:
                continue

            structure_builder.init_model(current_model_id, current_model_id)
            current_model_id += 1
            #this_model = model.get_id()
            #if current_model_id:
            #    structure_builder.init_model(current_model_id-1, this_model)
            #else:
            #    structure_builder.init_model(this_model)
            #prev_model = this_model

            for chain in model:
                if args.chain and chain != args.chain:
                    continue
                structure_builder.init_seg(' ')
                structure_builder.init_chain(chain.get_id())
                new_resseq = start_resseq
                for residue in chain:
                    hetero_flag, resseq, icode = residue.get_id()
                    if res_range_tuple and \
                            res_range_tuple[0] > resseq or \
                            res_range_tuple[1] <= resseq:
                        continue
                    new_resseq += 1
                    structure_builder.init_residue(residue.get_resname(), hetero_flag, new_resseq, icode)
                    for atom in residue:
                        structure_builder.init_atom(atom.get_name(), atom.get_coord(),
                                atom.get_bfactor(), atom.get_occupancy(), atom.get_altloc(),
                                atom.get_fullname())
                        structure_builder.set_line_counter(line_counter)
                        line_counter+= 1

    except FileNotFoundError:
        eprint("File not found: {}".format(structf))

    out_pdb = args.out
    out_structure = structure_builder.get_structure()
    pdbio = PDBIO(True)
    pdbio.set_structure(out_structure)
    pdbio.save(out_pdb, write_end=True)


if __name__ == "__main__":
    main()
