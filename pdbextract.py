#!/usr/bin/python3

from scan.biotools import *
import argparse

from Bio.PDB.StructureBuilder import StructureBuilder


def pdb_extract(structure, **kwargs):

    # model to extract from pdb
    extract_model = None if not 'model' in kwargs else kwargs['model']
    new_model_id = -1 if not 'new_model' in kwargs else kwargs['new_model']
    extract_chain = None if not 'chain' in kwargs else kwargs['chain']
    first_res = None if not 'first_res' in kwargs else kwargs['first_res']
    last_res = None if not 'last_res' in kwargs else kwargs['last_res']
    new_first_res = None if not 'new_first_res' in kwargs else kwargs['new_first_res']
    gap_count = 0 if not 'gap_count' in kwargs else kwargs['gap_count']
    water_id = None if not 'water' in kwargs else kwargs['water']

    model_rebumber_flag = bool( (extract_model is not None ) or new_model_id >= 0)
    res_renumber_flag = bool(first_res or last_res or new_first_res or gap_count)

    structure_builder = StructureBuilder()
    structure_builder.init_structure('pdb_extract')
    structure_builder.set_line_counter(0)
    line_counter = 0
    start_resseq_by_default = 1 if not new_first_res else new_first_res

    for model in structure:
        if model_rebumber_flag and \
                ( extract_model is not None ) and model.get_id() != extract_model:
            continue

        if model_rebumber_flag and new_model_id >= 0:
            this_model_id = new_model_id
            new_model_id += 1
        else:
            this_model_id = model.get_id()
        structure_builder.init_model(this_model_id, this_model_id)

        for chain in model:
            if extract_chain and chain.get_id() != extract_chain:
                continue
            structure_builder.init_seg(' ')
            structure_builder.init_chain(chain.get_id())
            resdict = {}
            if res_renumber_flag:
                # first_res = res_range_tuple[0]
                # last_res = res_range_tuple[1]

                resdict['before'] = select_residues_from_chain(chain,
                                        first_res = first_res, gap_count = gap_count)
                resdict['hit'] = select_residues_from_chain(chain,
                                        first_res = first_res, last_res = last_res)
                resdict['after'] = select_residues_from_chain(chain,
                                        last_res = last_res, gap_count = gap_count)
            else:
                resdict['before'] = []
                resdict['hit'] = chain.get_list()
                resdict['after'] = []
            new_resseq = start_resseq_by_default - len(resdict['before'])
            resdict['water'] = chain_water_id(chain, water_id)
            for key in ['before', 'hit', 'after', 'water']:
                for residue in resdict[key]:
                    if res_renumber_flag:
                        new_resid = ' ', new_resseq, ' '
                    else:
                        new_resid = residue.get_id()
                    structure_builder.init_residue(residue.get_resname(), *new_resid)
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
                    new_resseq += 1
                    if key == 'water' and gap_count and len(resdict['after']) != gap_count:
                        new_resseq += gap_count - len(resdict['after'])

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

    first_res, last_res  = None, None
    if args.residues:
        first_res, last_res = range_str_to_tuple(args.residues)

    try:
        structure = get_structure_from_file(args.struct, format='pdb')
        out_structure = pdb_extract(structure, model=args.model, chain=args.chain,
                    first_res=first_res, last_res=last_res, water=args.water)
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
