## DSSP
### Quick links:

- Home page: http://swift.cmbi.ru.nl/gv/dssp/
- Mac install: http://biskit.pasteur.fr/install/applications/dssp or http://proteinz.blogspot.com/2013/02/compiling-dssp-on-osx-lion-redux.html

### Installation for Ubuntu/Linux

1. Install the DSSP package

        sudo apt-get install dssp

2. The program installs itself as `mkdssp`, not `dssp`, and Biopython looks to execute `dssp`, so we need to symlink the name `dssp` to `mkdssp`
        
        sudo ln -s /usr/bin/mkdssp /usr/bin/dssp

3. Then you should be able to run `dssp`