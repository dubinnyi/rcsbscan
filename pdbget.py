#!/usr/bin/python3 -u

from scan.biotools import *
from numpy import pi
from Bio.PDB.vectors import Vector, calc_dihedral
import argparse


def main():
    arg_parser = argparse.ArgumentParser(
        usage='Extract data from set of structures')
    arg_parser.add_argument('struct', type = str, help = 'structures in pdb or mmcif format',
                            nargs = '+')
    arg_parser.add_argument('-r', '--recursive', action='store_true',
                            help='Recursive search of structures in folders', default=False)
    arg_parser.add_argument('-w', '--pdb-warnings', action='store_true',
                            help='show structure parsing warnings', default=False)
    arg_parser.add_argument('--phipsi', action='store_true',
                            help='Print phi/psi for residues', default=False)
    args = arg_parser.parse_args()

    if not args.pdb_warnings:
        warnings.simplefilter("ignore", PDBConstructionWarning)

    if args.recursive:
        struct_list = recursive_expand(args.struct)
        print("Recursive expansion: {} -> {} structures".format(len(args.struct),len(struct_list)))
    else:
        struct_list = args.struct

    nstruct = len(struct_list)
    print("nstruct: {}".format(nstruct))

    for structf in struct_list:
        try:
            (file_id, format) = parse_structure_filename(structf)
            structure = get_structure_from_file(structf)
            if not structure:
                continue
            for model in structure:
                for chain in model:
                    res_list = chain_sequence_aaselect(chain)
                    for r1,r2,r3 in zip(res_list[:-2],res_list[1:-1], res_list[2:]):
                        v1 = r1['C'].get_vector()
                        v2 = r2['N'].get_vector()
                        v3 = r2['CA'].get_vector()
                        v4 = r2['C'].get_vector()
                        v5 = r3['N'].get_vector()
                        phi = calc_dihedral(v1,v2,v3,v4)/pi*180
                        psi = calc_dihedral(v2,v3,v4,v5)/pi*180
                        print("{:6} {:3} {:1} {:4} {:4} PHI= {:8.3f} PSI= {:8.3f}".format(
                            structf,model.get_id(), chain.get_id(),
                            r2.get_id()[1], r2.get_resname(), phi,psi))

        except FileNotFoundError:
                eprint("File not found: {}".format(structf))


if __name__ == "__main__":
    main()
