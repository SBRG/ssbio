[![Binder](http://mybinder.org/badge.svg)](http://mybinder.org/repo/nmih/ssbio)

## ssbio
This Python package provides a collection of tools for people with questions in the realm
of structural systems biology. The main goals of this package are to:

1. Provide an easy way to map proteins to sequences and structures
2. Directly link structures to genome-scale SBML models
3. Prepare structures for downstream analyses, such as their use in molecular modeling software

##### Example questions you can answer with this package:

- How can I determine the number of protein structures available for my list of genes?
- What is the best, representative structure for my protein?
- Where, in a metabolic network, do these proteins work?
- Where do popular mutations show up on a protein?
- How can I compare the structural features of an entire proteome?
- and more...

## Quick install
Clone this repository to any location and then install it.

Cloning
```bash
$ git clone https://github.com/nmih/ssbio.git
```

Installation
```bash
$ pip uninstall ssbio  # Only if you have an old version installed
$ cd ssbio
$ python setup.py develop
```

Updating
```bash
$ cd ssbio
$ git pull
```

Uninstalling
```bash
pip uninstall ssbio
```

## Dependencies
ssbio heavily depends on Biopython, and for systems biology
applications COBRApy. If analyses are done in a Jupyter notebook,
Pandas DataFrames along with qgrid provide an easy way to look at
data. tqdm helps keep track of progress. Escher and NGLview allow for
interactive visualizations of metabolic maps and protein structures,
respectively.
