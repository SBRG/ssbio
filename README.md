[//]: [![Binder](http://mybinder.org/badge.svg)](http://mybinder.org/repo/nmih/ssbio)

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

## Documentation
Read the documentation [here](http://ssbio.readthedocs.io/en/latest/).

## Quick install

First install NGLview:
```bash
pip install nglview
```

Then install ssbio:
```bash
pip install ssbio
```

**Updating**
```bash
pip install ssbio --upgrade
```

**Uninstalling**
```bash
pip uninstall ssbio
```

### External programs to install
See: [Software Installations](https://github.com/SBRG/ssbio/wiki/Software-Installations)