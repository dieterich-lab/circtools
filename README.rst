**circtools**
======================================================================

**a one-stop software solution for circular RNA research**

.. figure:: https://raw.githubusercontent.com/dieterich-lab/circtools/master/docs/img/circtools_200px.png
   :alt: circtools

|docs| |build| |zenodo| |downloads| |pypi|

Introduction
-------------

Circular RNAs (circRNAs) originate through back-splicing events from linear primary transcripts, are resistant to exonucleases, typically not polyadenylated, and have been shown to be highly specific for cell type and developmental stage. Although few circular RNA molecules have been shown to exhibit miRNA sponge function, for the vast majority of circRNAs however, their function is yet to be determined.

The prediction of circular RNAs is a multi-stage bioinformatics process starting with raw sequencing data and usually ending with a list of potential circRNA candidates which, depending on tissue and condition may contain hundreds to thousands of potential circRNAs. While there already exist a number of tools for the prediction process (e.g. `DCC <https://github.com/dieterich-lab/DCC>`__ and `CircTest <https://github.com/dieterich-lab/CircTest>`__), publicly available downstream analysis tools are rare.

We developed **circtools**, a modular, Python3-based framework for circRNA-related tools that unifies several functionalities in single command line driven software. The command line follows the `circtools subcommand` standard that is employed in samtools or bedtools. Currently, circtools includes modules for detecting and reconstructing circRNAs,
a quick check of circRNA mapping results, RBP enrichment screenings, circRNA primer design, statistical testing, and an exon usage module.



Documentation
-------------

Click `here <http://docs.circ.tools/>`__ to access the complete documentation on Read the Docs.

Installation
------------

The ``circtools`` package is written in Python3 (>=3.4), two modules, namely ``detect`` and ``reconstruct`` also require a working Python 2 installation (>=2.7). It requires only a small number of external dependencies, namely standard bioinformatics tools:

-  `bedtools (>= 2.27.1) <https://bedtools.readthedocs.io/en/latest/content/installation.html>`__
   [RBP enrichment module, installed automatically]
-  `R (>= 3.3) <https://www.digitalocean.com/community/tutorials/how-to-install-r-on-ubuntu-16-04-2>`__
   [Data visualization and data processing]

Installation is managed through ``python3 setup.py install``. No sudo
access is required if the installation is executed with ``--user`` which
will install the package in a user-writeable folder. The binaries should
be installed to ``/home/$user/.local/bin/`` in case of Debian-based
systems.

``circtools`` was developed and tested on Debian Jessie but should also
run with any distribution.

The installation requires running python on the command line:

::

    git clone https://github.com/dieterich-lab/circtools.git
    cd circtools
    python3 setup.py install --verbose --user

The installation procedure will automatically install two dependencies:
`DCC <https://github.com/dieterich-lab/DCC>`__ and
`FUCHS <https://github.com/dieterich-lab/FUCHS>`__. The primer-design
module as well as the exon analysis and circRNA testing module require a
working installation of `R <https://cran.r-project.org/>`__ with
`BioConductor <https://www.bioconductor.org/install/>`__. All R packages
required are automatically installed during the setup. Please see the
`"Installing circtools" <http://docs.circ.tools/en/latest/Installation.html>`__
chapter of the main circtools documentation for more detailed installation instructions.

Modules
-------

Circtools currently offers seven modules:

detect `(detailed documentation) <https://circtools.readthedocs.io/en/latest/Detect.html>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``detect`` command is an interface to
`DCC <https://github.com/dieterich-lab/DCC>`__, also developed at the
Dieterich Lab. The module allows to detect circRNAs from RNA sequencing
data. The module is the foundation of all other steps for the circtools
work flow. All parameters supplied to circtools will be directly passed
to DCC.

quickcheck `(detailed documentation) <https://circtools.readthedocs.io/en/latest/Quickcheck.html>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The quickcheck module of circtools is an easy way to check the results
of a DCC run for problems and to quickly assess the number of circRNAs
in a given experiment. The module needs the mapping log files produced
by STAR as well as the directory with the DCC results. The module than
generates a series of figures in PDF format to assess the results.

reconstruct `(detailed documentation) <https://circtools.readthedocs.io/en/latest/Reconstruct.html>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``reconstruct`` command is an interface to
`FUCHS <https://github.com/dieterich-lab/FUCHS>`__. FUCHS is employing
DCC-generated data to reconstruct circRNA structures. All parameters
supplied to circtools will be directly passed to FUCHS.

circtest `(detailed documentation) <https://circtools.readthedocs.io/en/latest/Circtest.html>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``circtest`` command is an interface to
`CircTest <https://github.com/dieterich-lab/CircTest>`__. The module a a
very convenient way to employ statistical testing to circRNA candidates
generated with DCC without having to write an R script for each new
experiment. For detailed information on the implementation itself take a
look at the `CircTest
documentation <https://github.com/dieterich-lab/CircTest>`__. In
essence, the module allows dynamic grouping of the columns (samples) in
the DCC data.

exon `(detailed documentation) <https://circtools.readthedocs.io/en/latest/Exon.html>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The exon module of circtools employs the `ballgown R
package <https://www.bioconductor.org/packages/release/bioc/html/ballgown.html>`__
to combine data generated with DCC and circtest with ballgown-compatible
``stringtie`` output or cufflinks output converted via
`tablemaker <https://github.com/leekgroup/tablemaker>`__ in order get
deeper insights into differential exon usage within circRNA candidates.

enrich `(detailed documentation) <https://circtools.readthedocs.io/en/latest/Enrichment.html>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``enrichment`` module may be used to identify circRNAs enriched for
specific RNA binding proteins (RBP) based on DCC-identified circRNAs and
processed
`eCLIP <http://www.nature.com/nmeth/journal/v13/n6/full/nmeth.3810.html>`__
data. For K526 and HepG2 cell lines plenty of this data is available
through the
`ENCODE <https://www.encodeproject.org/search/?type=Experiment&assay_title=eCLIP>`__
project.

primer `(detailed documentation) <https://circtools.readthedocs.io/en/latest/primer.html>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``primer`` command is used to design and visualize primers required
for follow up wet lab experiments to verify circRNA candidates.

.. |docs| image:: https://readthedocs.org/projects/circtools/badge/?version=latest
    :alt: Documentation Status
    :scale: 100%
    :target: https://circtools.readthedocs.io/en/latest/?badge=latest

.. |build| image:: https://travis-ci.org/dieterich-lab/circtools.svg?branch=master
    :alt: Build Status
    :scale: 100%
    :target: https://travis-ci.org/dieterich-lab/circtools

.. |zenodo| image:: https://zenodo.org/badge/83248654.svg
    :alt: Zenodo DOI link
    :scale: 100%
    :target: https://zenodo.org/badge/latestdoi/83248654

.. |downloads| image:: https://pepy.tech/badge/circtools
    :alt: Python Package Index Downloads
    :scale: 100%
    :target: https://pepy.tech/project/circtools

.. |pypi| image:: https://badge.fury.io/py/circtools.svg
    :alt: Python package version
    :scale: 100%
    :target: https://badge.fury.io/py/circtools
