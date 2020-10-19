#!/usr/bin/python3 -u

from scan.biotools import *
from scan.struct4fit import Struct4Fit
from scan.fitcounter import FitCounter

import argparse


def pool_fit_one(tuple_of_args):
    (ref_struct, one_file, counter) = tuple_of_args
    try:
        counter.timer_start()
        ref_struct.scan(one_file, counter)
    except FileNotFoundError:
        eprint("Failed to scan structure {}".format(one_file))
        counter.new_error()
    counter.timer_finish()
    return counter


def main():
    arg_parser = argparse.ArgumentParser(
        usage='Fit reference structure to every substructute')
    arg_parser.add_argument('struct', type = str, help = 'structures in pdb or mmcif format', nargs = '+')
    arg_parser.add_argument('-p', '--print-header', action='store_true',
                            help='print header')
    arg_parser.add_argument('--ref-structure', type=str, help='Reference structure', required=True)
    arg_parser.add_argument('--ref-model', type=int, help='Model number in reference structure', default = None)
    arg_parser.add_argument('--ref-chain', type=str, help='Chain in reference structure', default = None)
    arg_parser.add_argument('--ref-residues', type=str, help='Residue range in reference structure', default = None)
    arg_parser.add_argument('--ref-atoms', type=str, help='Atoms in reference structure',
                            default='N,CA,C,O')
    arg_parser.add_argument('-w', '--pdb-warnings', action='store_true',
                            help='show structure parsing warnings', default=False)
    arg_parser.add_argument('-v', '--verbose', action='store_true',
                            help='Verbose output', default=False)
    arg_parser.add_argument('-r', '--recursive', action='store_true',
                            help='Recursive search of structures in folders', default=False)
    arg_parser.add_argument('--max-rms', type=float, help='Maximum RMSD to print, ingore otherwise', default=1.0)
    arg_parser.add_argument('--water', type=str, help='Water molecule to scan (residue number)')
    arg_parser.add_argument('--water-max-rms', type=float, help='Max rms for water match', default=2.0)
    arg_parser.add_argument('--save-pdb-hits', type=str, help='Save pdb hits to file', default = None)
    arg_parser.add_argument('--renumber-pdb', action='store_true',
                            help='Renumber residues in pdb hits', default=False)
    arg_parser.add_argument('--xray-res', type = float,
                            help='Maximal resolution of X-ray structures to scan')
    arg_parser.add_argument('--xray-only', action='store_true', default=False,
                            help='Scan only X-ray structures')
    arg_parser.add_argument('--ncpu', type=int, help='Number of CPU to use', default = None)
    args = arg_parser.parse_args()

    if not args.pdb_warnings:
        warnings.simplefilter("ignore", PDBConstructionWarning)

    #ref_struct=get_structure_from_file(args.ref_structure)

    s4fit = Struct4Fit(args.ref_structure, args.ref_model, args.ref_chain,
                       args.ref_residues,  args.ref_atoms, args.verbose, args.max_rms)

    if args.water:
        s4fit.set_water(args.water, args.water_max_rms)

    if not s4fit:
        quit(-1)

    print("save pdb hits: {}".format(args.save_pdb_hits))
    if args.save_pdb_hits:
        s4fit.set_write(args.save_pdb_hits)

    if args.xray_res:
        s4fit.set_Xray(args.xray_res, args.xray_only)

    print(str(s4fit))

    ##############################
    #  Start fit scan
    ##############################

    print("Start fit scan")

    CountTotal = FitCounter(s4fit.result_name)
    if args.verbose:
        print(CountTotal)

    ncpu_onboard = mp.cpu_count()
    if args.ncpu:
        ncpu = min(ncpu_onboard, args.ncpu)
    else:
        ncpu = ncpu_onboard
    print("Start the pool of {} CPU (of {} available)".format(ncpu, ncpu_onboard))
    pool = mp.Pool(ncpu)

    if not args.struct:
        # Should not be here: nargs='+'
        eprint("No structures selected. Exiting")
        quit(-1)

    print("struct list: {}".format(args.struct))

    if args.recursive:
        struct_list = recursive_expand(args.struct)
        print("Recursive expansion: {} -> {} structures".format(len(args.struct),len(struct_list)))
    else:
        struct_list = args.struct

    nstruct = len(struct_list)
    print("nstruct: {}".format(nstruct))

    print("Prepare arguments for {} structures".format(nstruct))
    all_args = [(s4fit, struct_list[i], FitCounter(s4fit.result_name))
                for i in range(nstruct)]

    print("Start map_async".format(nstruct))
    walltime_start = time.monotonic()
    async_result = pool.map_async(pool_fit_one, all_args)

    print("map_async submitted {} tasks to Pool of {} cpu".format(nstruct, ncpu))
    pool.close()
    pool.join()
    async_result.wait()
    walltime_finish = time.monotonic()
    print("Pool finished evaluation of {} tasks".format(nstruct))
    result_counters = async_result.get()
    # print("Results obtained".format(nstruct))
    #for structf in args.struct:
    #    try:
    #        s4fit.scan(structf, c)
    #    except FileNotFoundError:
    #        eprint("Failed to read structure {}".format(structf))
    #        continue

    for rc in result_counters:
        if args.verbose:
            print(rc)
        CountTotal = CountTotal + rc

    walltime = walltime_finish - walltime_start
    print("All counters are accumulated")
    CountTotal.sort_hits_by_rmsd()
    print("All hits are sorted by RMSD")
    CountTotal.print_hits()
    print("Overall statistics:")
    print(CountTotal)
    print("Evaluation time: {:8.2f}".format(walltime))
    print("Save all hits to {}".format(args.save_pdb_hits))
    CountTotal.save_hits_to_pdb(args.save_pdb_hits)



if __name__ == "__main__":
    main()
