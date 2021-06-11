**rcsbscan** is a biopython-based tool for large-scale RCSB protein 
structure database screening with a user-defined 3D structural template (reference structure).

The amino acid sequence of the template is completely ignored in the scan process, only atomwise RMSD of the reference 
structure with residue tuples is the criteria for the structural match hit. Residue tuple is the sequance of
residues with defined 3D structure. All rmsd hits are prinded with match description and
are optionally rwitten to the separate PDB file for future analysis.

The RCSB protein databease should be downloaded locally and preferably to the hard drive of the 
computational server. Please follow the instructions provided by RCSB on the following link:
https://www.rcsb.org/docs/programmatic-access/batch-downloads-with-shell-script


```
rcsbscan.py ARGS 
    Scan database by sequential atomwise fit of the provided  reference 
    structure to every substructute

Positional arguments:
  struct                Structure(s) in pdb, mmcif or mmtf format, optionally gzipped. 
                        Directory or a list of directories with protein structures is allowed with -r flag (see below).
                        To scan full RCSB database, first download it locally and privide it with an additional -r option

Optional arguments:
  -h, --help            show this help message and exit
  -p, --print-header    Print header of the PDB which is scanned
  --ref-structure REF_STRUCTURE
                        Reference structure, in PDB, MMCIF of MMTF formats
  --ref-model REF_MODEL
                        Model number in the reference structure, e.g. '0'
  --ref-chain REF_CHAIN
                        Chain in the reference structure, e.g. 'A'
  --ref-residues REF_RESIDUES
                        Residue range in the reference structure, e.g. '26-33'
  --ref-atoms REF_ATOMS
                        Atoms in reference structure. 
                        By default, only four atoms per residue are considured: N, CA, C, O
  -w, --pdb-warnings    Show structure parsing warnings
  -v, --verbose         Verbose output
  -r, --recursive       Recursive search of structures in folders
  --max-rms MAX_RMS     Maximum RMSD to print [ default 1.0 A ] 
  --water WATER         Water molecule to be included to scan (residue number)
  --water-max-rms WATER_MAX_RMS
                        Max rms for water match [ default 2.0 A ] 
  --save-pdb-hits SAVE_PDB_HITS
                        Save pdb hits to file
  --renumber-pdb        Renumber residues in the output pdb hits
  --xray-res XRAY_RES   Maximal resolution of X-ray structures to scan
  --xray-only           Scan only X-ray structures
  --ncpu NCPU           Number of CPU to use [ All available CPU's are used by default ] 
```
**Examples:**

Ussuming that RCSB clone is downloadad to the folder ```./RCSB/pdb/```:
```
rcsbscan.py ./RCSB/pdb/ -r --ref-structure ./examples/alpha-helix_A10.pdp
```
Scans every PDB files of the local RCSB clone for atomwise match with an alpha-helix template
