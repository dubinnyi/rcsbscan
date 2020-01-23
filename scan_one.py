#!/usr/bin/python

import argparse
import os
import gzip


def get_pdb_sequence_info(filename):
   if filename[-2:] == "gz":
      openfunc = gzip.open
      openarg = 'rt'
   else:
      openfunc = open
      openarg = 'r'
   seq_info = {}
   with openfunc(filename, openarg) as f:
      lines = f.readlines()
      for l in lines:
         if l.startswith("SEQRES   1 "):
            seq_info['length'] = int(l[13:18])
         elif l.startswith("ATOM ") or l.startswith("HETATM "):
            return(seq_info)
   return seq_info

def get_pdb_sequence_length(filename):
   seq_info = get_pdb_sequence_info(filename)
   return(seq_info['length'])

def main():
   parser = argparse.ArgumentParser(
      usage='scan_one [options] pdbfile.gz')
   parser.add_argument('pdb', type=str, help='structures in .pdb.gz format')
   parser.add_argument('--max-rmsd', type=float, help='max rmsd to print', default=1.5)
   args = parser.parse_args()

   sub_length =8
   file_to_scan = args.pdb
   cmd_base="./calculate_rmsd.py 4OM4_A.pdb.gz " + \
      " {pdb} --select --res-range-b {{i}}-{{j}} --atom-names N,CA,C,O " + \
      " --max-rmsd {rmsd} "
   cmd_base = cmd_base.format(pdb = file_to_scan, rmsd = args.max_rmsd)
   
   seq_length = get_pdb_sequence_length(file_to_scan)
   
   print("{} seq_length={}".format(file_to_scan, seq_length))
   print(cmd_base)

   for i in range(1,seq_length+1-(sub_length-1)):
      j=i+sub_length-1
      cmd=cmd_base.format(i = i, j = j)
      #print(cmd)
      os.system(cmd)
  
if __name__ == "__main__":
    main()
