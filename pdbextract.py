#!/usr/bin/python3

from scan.biotools import *
import argparse

from Bio.PDB.StructureBuilder import StructureBuilder


def pdb_extract(structure, **kwargs):

    extract_model = None if not 'model' in kwargs else kwargs['model']
    current_model_id = 0 if not 'new_model' in kwargs else kwargs['new_model']
    extract_chain = None if not 'chain' in kwargs else kwargs['chain']
    res_range_tuple = None if not 'res_range' in kwargs else kwargs['res_range']
    water_id = None if not 'water' in kwargs else kwargs['water']

    structure_builder = StructureBuilder()
    structure_builder.init_structure('pdb_extract')
    structure_builder.set_line_counter(0)
    line_counter = 0
    start_resseq = -1

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
            first_res = res_range_tuple[0]
            last_res = res_range_tuple[1]
            resdict = {}
            resdict['before'] = select_residues_from_chain(chain,
                                                           first = first_res, count = 1)
            resdict['hit'] = select_residues_from_chain(chain,
                                                        first=first_res, last=last_res)
            resdict['after'] = select_residues_from_chain(chain,
                                                          last = last_res, count = 1)
            resdict['water'] = chain_water_id(chain, water_id)
            for key in ['before', 'hit', 'after', 'water']:
                for residue in resdict[key]:
                    new_resseq += 1
                    structure_builder.init_residue(residue.get_resname(), ' ', new_resseq, ' ')
                    residue_atoms = None
                    if key == 'before':
                        residue_atoms = [atom for atom in residue if \
                                         (atom.get_name() == 'C' or atom.get_name() == 'O')]
                    elif key == 'hit' or key == 'water':
                        residue_atoms = residue.get_list()
                    elif key == 'after':
                        residue_atoms = [atom for atom in residue if \
                                         (atom.get_name() == 'N' or atom.get_name() == 'HN')]
                    for atom in residue_atoms:
                        structure_builder.init_atom(atom.get_name(), atom.get_coord(),
                                                    atom.get_bfactor(), atom.get_occupancy(),
                                                    atom.get_altloc(),
                                                    atom.get_fullname())
                        structure_builder.set_line_counter(line_counter)
                        line_counter += 1

            # for residue in res_before:
            #     residue_atoms = [atom for atom in residue if \
            #                      (atom.get_name() == 'C' or atom.get_name() == 'O')]
            #     new_resseq += 1
            #     structure_builder.init_residue(residue.get_resname(), ' ', new_resseq, ' ')
            #     for atom in residue_atoms:
            #         structure_builder.init_atom(atom.get_name(), atom.get_coord(),
            #                                     atom.get_bfactor(), atom.get_occupancy(),
            #                                     atom.get_altloc(),
            #                                     atom.get_fullname())
            #         structure_builder.set_line_counter(line_counter)
            #         line_counter += 1
            #
            # for residue in chain:
            #     residue_atoms = None
            #     resid = residue.get_id()
            #     hetero_flag, resseq, icode = resid
            #     if res_range_tuple and ( not hetero_flag or hetero_flag == ' ' ) :
            #         if res_range_tuple[0] - 1 > resseq:
            #             continue
            #         elif res_range_tuple[0] - 1 == resseq:
            #             residue_atoms = [ atom for atom in residue if \
            #                              ( atom.get_name() == 'C' or atom.get_name() == 'O' ) ]
            #         elif res_range_tuple[1] < resseq:
            #             continue
            #         elif res_range_tuple[1] == resseq:
            #             residue_atoms = [ atom for atom in residue if \
            #                              (atom.get_name() == 'N' or atom.get_name() == 'HN' ) ]
            #     if not residue_atoms:
            #         residue_atoms = residue.get_list()
            #     if hetero_flag == 'W':
            #         if not water_id:
            #             continue
            #         elif resseq != water_id:
            #             continue
            #     if hetero_flag.startswith('H'):
            #         continue
            #     new_resseq += 1
            #     structure_builder.init_residue(residue.get_resname(), hetero_flag, new_resseq, icode)
            #     for atom in residue_atoms:
            #         structure_builder.init_atom(atom.get_name(), atom.get_coord(),
            #                                     atom.get_bfactor(), atom.get_occupancy(), atom.get_altloc(),
            #                                     atom.get_fullname())
            #         structure_builder.set_line_counter(line_counter)
            #         line_counter += 1
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
        out_structure = pdb_extract(structure, model=args.model, chain=args.chain,
                            res_range=res_range_tuple, water=args.water)
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
