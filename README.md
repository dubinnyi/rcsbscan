**rcsbscan** is a biopython-based tool for large-scale RCSB protein 
structure database screening with a user-defined 3D structural template (reference structure).

The amino acid sequence of the template is completely ignored in the scan process, only atomwise RMSD of the reference 
structure with **residue tuples** by all backbone atoms is the criteria for the structural hit. 
**Residue tuple** is the sequence of residues with defined 3D structure of the same length as the 
reference structural template provided. For example, 9-residue protein have four six-reside tuples to be checked for 
structural hit:
```buildoutcfg
1-ABCDEFGHI-9  -- scanned structure from the database
  XXXXXX  |
  |XXXXXX |
  | XXXXXX|
  |  XXXXXX    -- reference structure
 
```
All rmsd hits are printed to the terminal with match description. The structures found 
are optionally saved to the separate PDB file for future analysis.

The RCSB protein databease should be downloaded locally and preferably to the hard drive of the 
computational server. Please follow the instructions provided by RCSB on the following link:
https://www.rcsb.org/docs/programmatic-access/batch-downloads-with-shell-script

```
rcsbscan.py ARGS 
    Scan database by sequential atomwise fit of the provided  reference 
    structure to every substructute

Positional arguments:
  struct                Structure(s) in pdb, mmcif or mmtf format, optionally gzipped. 
                        Directory or a list of directories with protein structures is 
                        allowed with -r flag (see below). To scan full RCSB database, 
                        first download it locally and privide it with an additional 
                        -r option

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
                        By default, only four atoms per residue are considured: 
                        N, CA, C, O
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

Ussuming that RCSB clone is downloadad to the folder ```~/RCSB/pdb/```:
```
rcsbscan.py ~/RCSB/pdb/m0/pdb1m0k.ent.gz --ref-structure ./examples/alpha-helix_A10.pdp
```
-- scans the ``1M0K`` structure of bacteriorhodopsin at 1.43 A resolution 
( (2002) J Mol Biol 321: 715-726 )
atomwise and prints all structural hits with an alpha-helix template provided in the example file. 
The output will find numerous alpha helices in the structure of bacteriorhodopsin:
```buildoutcfg
save pdb hits: None
REF_4FIT: ALPHA-HELIX_A10 model=   0, chain=   seq=    1 AAAAAAAAAA 10   atoms=N,CA,C,O fit_atoms=40 max_rms=1.0000
Start fit scan
Start the pool of 8 CPU (of 8 available)
struct list: ['/home/maxim/RCSB/pdb/m0/pdb1m0k.ent.gz'] (trancated at 10 structures)
nstruct: 1
Prepare arguments for 1 structures
Start map_async
map_async submitted 1 tasks to Pool of 8 cpu
RMSD_HIT: 1M0K XRay 1.43A model=   0, chain= A size=  222 hit=    8 PEWIWLALGT 17   rms= 0.8595
RMSD_HIT: 1M0K XRay 1.43A model=   0, chain= A size=  222 hit=    9 EWIWLALGTA 18   rms= 0.3121
RMSD_HIT: 1M0K XRay 1.43A model=   0, chain= A size=  222 hit=   10 WIWLALGTAL 19   rms= 0.2138
RMSD_HIT: 1M0K XRay 1.43A model=   0, chain= A size=  222 hit=   11 IWLALGTALM 20   rms= 0.2530
RMSD_HIT: 1M0K XRay 1.43A model=   0, chain= A size=  222 hit=   12 WLALGTALMG 21   rms= 0.2440
....
1M0K XRay 1.43A model=   1, chain= A size=  222 hit=   81 ARYADWLFTT 90   rms= 0.8996
1M0K XRay 1.43A model=   0, chain= A size=  222 hit=   82 RYADWLFTTP 91   rms= 0.9623
1M0K XRay 1.43A model=   1, chain= A size=  222 hit=   82 RYADWLFTTP 91   rms= 0.9623
Overall statistics:
fitscan statistics after    3 sec:
        1 files with structutes
        0 files skipped
        1 files scanned
        2 structures in all models/chains
      228 hits
      406 tuples of 10 residues superimposed and rms of atoms N,CA,C,O evaluated
        0 errors
Evaluation time:     3.82
No filename was provided to store the results found
Use the command-line argimens '--save-pdb-hits FILE' to save all matches in PDB format
```
