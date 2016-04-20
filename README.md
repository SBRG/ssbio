## Install
```bash
$ python setup.py install
```


## Clean PDB files
### Description
This script will automatically:
    * Add missing chains to a PDB file
    * Select a single chain if noted
    * Remove alternate atom locations
    * Add atom occupancies
    * Add B (temperature) factors (default Biopython behavior)
Cleaned PDBs will be in a clean_pdbs folder where the script is executed.
#### Example: script help
```bash
$ cleanpdb.py --help
```
#### Example: clean one PDB file
```bash
$ cleanpdb.py 1kf6.pdb
```
#### Example: clean multiple PDB files
```bash
$ cleanpdb.py *.pdb
```
#### Example: clean a whole directory of PDB
```bash
$ cleanpdb.py /path/to/pdb/files
```

## Mutate a PDB file
### Description
This script will automatically:
    * Mutate one PDB file with a specified chain and residue number to a new amino acid
    * Remove the heavy atoms of the old residue
The mutated PDB will be in a mutated_pdbs folder. The residue will only change names to the mutated residue, heavy atoms will not be added (run tleap.py for that!).

Specify your mutation in the form of:
```
CHAIN.RESNUM.MUTATEAA [, CHAIN.RESNUM.MUTATEAA, ...]
```
Where CHAIN is the chain ID in which RESNUM is contained, which will be mutated to MUTATEAA. MUTATEAA can be in the form of a 1 or 3 letter code (ie. T or Tyr).
#### Example: script help
```bash
$ mutatepdb.py --help
```
#### Example
This will mutate residue 3 in chain A to an alanine, and residue 3 in chain B to an alanine.
```bash
$ mutatepdb.py 1kf6.pdb A.3.ALA,B.3.ALA
```
#### Example with cleaning option
This will keep only the chains which are being mutated, as well as do the default cleaning options (required for tleap to run properly)
```bash
$ mutatepdb.py 1kf6.pdb A.3.ALA,B.3.ALA --clean
```

## Run tleap to prepare files for AMBER MD
### Description
This script will automatically:
    * Run tleap for all .pdb files in a folder
    * Create a xleap_modified and amber_minimized folder with files ready for MD simulation
Caution: this will most likely fail if PDB files are not cleaned.
### Requirements
AMBERTools needs to be installed on your computer for tleap to run properly. For instructions see: http://jswails.wikidot.com/installing-amber14-and-ambertools14
Test the installation by running
```
$ tleap
```
from your terminal to make sure that AMBERTools properly.
#### Example: script help
```bash
$ tleap.py --help
```
#### Example
This will run tleap on all files in a folder
```bash
$ tleap.py .
```