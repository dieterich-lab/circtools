# **circtools** - a one-stop software solution for circular RNA research

![circtools](docs/img/circtools_500px.png)

## Documentation

Read the Docs: [![Documentation Status](https://readthedocs.org/projects/circtools/badge/?version=latest)](http://circtools.readthedocs.io/en/latest/?badge=latest) or just click [here](http://circtools.readthedocs.io/en/latest/) to access the complete documentation.

## Installation

This package is written in python3 (3.4). It only a small number of external dependencies, namely standard bioinformatics tools:

* [bedtools (>= 2.27.1)](http://bedtools.readthedocs.io/en/latest/content/installation.html) [RBP enrichment module, installed automatically]
* [R (>= 3.3)](https://www.digitalocean.com/community/tutorials/how-to-install-r-on-ubuntu-16-04-2) [Data visualization and data processing] 

Installation is managed through `python3 setup.py install`. No sudo access is required if the installation is executed with ``--user`` which will install the package in a user-writeable folder. The binaries should be installed to ``/home/$user/.local/bin/`` in case of Debian-based systems.

``circtools`` was developed and tested on Debian Jessie but should also run with any distribution.

The installation requires running python on the command line:

```
git clone git@github.com:dieterich-lab/circtools.git
cd circtools
python3 setup.py install --verbose --user
```

The installation procedure will automatically install two dependencies: [DCC](https://github.com/dieterich-lab/DCC) and [FUCHS](https://github.com/dieterich-lab/FUCHS). The primer-design module as well as the exon analysis and circRNA testing module require a working installation of [R](https://cran.r-project.org/) with [BioConductor](https://www.bioconductor.org/install/). All R packages required are automatically installed during the setup. Please see [Installing circtools](http://circtools.readthedocs.io/en/latest/Installation.html) for more detailed installation instructions. A more detailed introduction to circtools can be found [on the readthedocs pages](http://circtools.readthedocs.io/en/latest/index.html) 

## Modules

Circtools currently offers seven modules:

### detect [(detailed documentation)](http://circtools.readthedocs.io/en/latest/Detect.html)

The ``detect`` command is an interface to [DCC](https://github.com/dieterich-lab/DCC), also developed at the Dieterich Lab. The module allows to detect circRNAs from RNA sequencing data. The module is the foundation of all other steps for the circtools work flow. All parameters supplied to circtools will be directly passed to DCC.

### quickcheck [(detailed documentation)](http://circtools.readthedocs.io/en/latest/Quickcheck.html)

The quickcheck module of circtools is an easy way to check the results of a DCC run for problems and to quickly assess the number of circRNAs in a given experiment. The module needs the mapping log files produced by STAR as well as the directory with the DCC results. The module than generates a series of figures in PDF format to assess the results.

### reconstruct [(detailed documentation)](http://circtools.readthedocs.io/en/latest/Reconstruct.html)

The ``reconstruct`` command is an interface to [FUCHS](https://github.com/dieterich-lab/FUCHS). FUCHS is employing DCC-generated data to reconstruct circRNA structures. All parameters supplied to circtools will be directly passed to FUCHS.

### circtest [(detailed documentation)](http://circtools.readthedocs.io/en/latest/Circtest.html)

The ``circtest`` command is an interface to [CircTest](https://github.com/dieterich-lab/CircTest). The module a a very convenient way to employ statistical testing to circRNA candidates generated with DCC without having to write an R script for each new experiment. For detailed information on the implementation itself take a look at the [CircTest documentation](https://github.com/dieterich-lab/CircTest). In essence, the module allows dynamic grouping of the columns (samples) in the DCC data. 

### exon [(detailed documentation)](http://circtools.readthedocs.io/en/latest/Exon.html)
 
The exon module of circtools employs the [ballgown R package](https://www.bioconductor.org/packages/release/bioc/html/ballgown.html) to combine data generated with DCC and circtest with ballgown-compatible `stringtie` output or cufflinks output converted via [tablemaker](https://github.com/leekgroup/tablemaker) in order get deeper insights into differential exon usage within circRNA candidates. 

### enrich [(detailed documentation)](http://circtools.readthedocs.io/en/latest/Enrichment.html)

The ``enrichment`` module may be used to identify circRNAs enriched for specific RNA binding proteins (RBP) based on DCC-identified circRNAs and processed [eCLIP](http://www.nature.com/nmeth/journal/v13/n6/full/nmeth.3810.html) data. For K526 and HepG2 cell lines plenty of this data is available through the [ENCODE](https://www.encodeproject.org/search/?type=Experiment&assay_title=eCLIP)
 project. The enrich module understands the following options:
 
### primer [(detailed documentation)](http://circtools.readthedocs.io/en/latest/primer.html)

The ``primer`` command is used to design and visualize primers required for follow up wet lab experiments to verify circRNA candidates.




