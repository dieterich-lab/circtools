This project contains the framework of the circular RNA toolbox ``circtools``.

# Installation

This package is written in python3 (3.4). It has a number of external dependencies, mostly standard bioinformatics tools:

* [bedtools (>= 2.26.0)](http://bedtools.readthedocs.io/en/latest/content/installation.html) [RBP enrichment module, installed automatically]
* [R (>= 3.3)](https://www.digitalocean.com/community/tutorials/how-to-install-r-on-ubuntu-16-04-2)
  [Primer design module]
* [OligoArrayAux](http://unafold.rna.albany.edu/?q=DINAMelt/OligoArrayAux)
  [required by DECIPHER Bioconductor package for annealing efficiency estimations, installed automatically]

Installation is managed through `python3 setup.py install`. No sudo access is required if the installation is executed with ``--user`` which will install the package in a user-writeable folder. The binaries should be installed to ``/home/$user/.local/bin/`` in case of Debian-based systems.

``Circtools`` was developed and tested on Debian Jessie.

The installation requires running python on the command line:

```
git clone git@github.com:dieterich-lab/circtools.git
cd circtools
python3 setup.py install --verbose --user
```

The installation procedure will automatically install two dependencies: [DCC](https://github.com/dieterich-lab/DCC) and [FUCHS](https://github.com/dieterich-lab/FUCHS). The primer-design module as well as the exon analysis and circRNA testing module require a working installation of [R](https://cran.r-project.org/) with [BioConductor](https://www.bioconductor.org/install/). All R packages required are automatically installed during the setup.

