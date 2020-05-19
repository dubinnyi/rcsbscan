#!/usr/bin/python3 -u

from scan.biotools import *
import numpy as np
from Bio.PDB.vectors import Vector, calc_dihedral
import argparse


def main():
    arg_parser = argparse.ArgumentParser(
        usage='Extract data from set of structures')
    arg_parser.add_argument('struct', type = str, help = 'structures in pdb or mmcif format',
                            nargs = '+')
    arg_parser.add_argument('-o', '--output', type=str, help='Output file name')
    arg_parser.add_argument('-r', '--recursive', action='store_true',
                            help='Recursive search of structures in folders', default=False)
    arg_parser.add_argument('-w', '--pdb-warnings', action='store_true',
                            help='show structure parsing warnings', default=False)
    arg_parser.add_argument('--phipsi', action='store_true',
                            help='Print phi/psi for residues', default=False)
    arg_parser.add_argument('-v', '--verbose', action='store_true',
                            help='Verbose output', default=False)
    args = arg_parser.parse_args()

    if not args.pdb_warnings:
        warnings.simplefilter("ignore", PDBConstructionWarning)

    if args.recursive:
        struct_list = recursive_expand(args.struct)
        print("Recursive expansion: {} -> {} structures".format(len(args.struct),len(struct_list)))
    else:
        struct_list = args.struct

    nfiles = len(struct_list)
    print("nfiles: {}".format(nfiles))

    phi_array, psi_array, info = [], [], []
    maxres = 0
    for structf in struct_list:
        try:
            structure = get_structure_from_file(structf)
            if not structure:
                eprint("No structure found in file {}".format(structf))
                continue
            print("Read structure(s) from file {}".format(structf))
            s = 0
            for model in structure:
                for chain in model:
                    res_list = chain_sequence_aaselect(chain)
                    Nres = len(res_list) - 2
                    if Nres > maxres:
                        maxres = Nres
                        print("maxres = {} in file {}, structure {}".format(maxres, structf, s))
                    #maxres = Nres if Nres > maxres else maxres
                    phi_list, psi_list = [], []
                    for res1,res2,res3 in zip(res_list[:-2],res_list[1:-1], res_list[2:]):
                        v1 = res1['C'].get_vector()
                        v2 = res2['N'].get_vector()
                        v3 = res2['CA'].get_vector()
                        v4 = res2['C'].get_vector()
                        v5 = res3['N'].get_vector()
                        phi = calc_dihedral(v1, v2, v3, v4)/np.pi*180
                        psi = calc_dihedral(v2, v3, v4, v5)/np.pi*180
                        phi_list.append(phi)
                        psi_list.append(psi)
                        info_str = "{:6} {:3} {:1} {:4} {:4} PHI= {:8.3f} PSI= {:8.3f}".format(
                                structf,model.get_id(), chain.get_id(),
                                res2.get_resname(), res2.get_id()[1] + res2.get_id()[2],
                                phi, psi)
                        info.append(info_str)
                        if args.verbose:
                            print(info_str)
                    phi_array.append(phi_list)
                    psi_array.append(psi_list)
                    s += 1

        except FileNotFoundError:
                eprint("File not found: {}".format(structf))

    Nstruct = len(phi_array)
    print("All {} structures read, maxres = {}".format(Nstruct, maxres))
    phi_np = np.zeros((Nstruct, maxres))
    psi_np = np.zeros((Nstruct, maxres))
    #phi_np = np.array(phi_array)
    #psi_np = np.array(psi_array)
    for s in range(Nstruct):
        try:
            phi_np[s,:] = phi_array[s]
            psi_np[s,:] = psi_array[s]
        except:
            eprint("Can't copy phi,psi from structure {}:".format(s))
            eprint(info[s])
            continue
    print("PHI shape = {}".format(phi_np.shape))
    print("PSI shape = {}".format(psi_np.shape))
    if args.output:
        np.save(args.output + '_phi.npy', phi_np)
        np.save(args.output + '_psi.npy', psi_np)



if __name__ == "__main__":
    main()
