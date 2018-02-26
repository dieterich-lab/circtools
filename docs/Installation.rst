Installation
=============


**circtools** is written in Python3 (>=3.4). The tool has a number of external dependencies, mostly standard bioinformatics tools and packages. The installation will, by default, try to install all required dependencies.

Installation is performed via `python3 setup.py install`. No sudo access is required if the installation is suffixed with ``--user`` which will install the package in a user-writeable folder. In this case, the binaries should be installed to ``/home/$user/.local/bin/`` (for Debian-based systems).

``circtools`` was developed and tested on Debian Jessie 8 64 Bit.


Quick 3-line installation
--------------------------

The default installation requires running python on the command line and will install everything needed to run circtools *except bedtools and R* (see below)::

    git clone git@github.com:dieterich-lab/circtools.git
    cd circtools
    python3 setup.py install --verbose --user


Required dependencies
---------------------

External tools
^^^^^^^^^^^^^^^


* `bedtools [>= 2.27.0] <http://bedtools.readthedocs.io/en/latest/content/installation.html>`_ required by the enrichment module

.. important:: Currently at least bedtools v2.26.0-125-g52db654-dirty is required to run the enrichment module. Earlier versions have issues with the ``shuffle`` command.



* `R [>= 3.3] <https://www.digitalocean.com/community/tutorials/how-to-install-r-on-ubuntu-16-04-2>`_ required by visualisation scripts and the primer design module

* `OligoArrayAux <http://unafold.rna.albany.edu/?q=DINAMelt/OligoArrayAux>`_ required by DECIPHER Bioconductor package for annealing efficiency estimations, *installed automatically*


The installation procedure will automatically install two additional Python-based dependencies: `DCC <https://github.com/dieterich-lab/DCC>`_ and `FUCHS <https://github.com/dieterich-lab/FUCHS>`_ by temporarily cloning the repositories and installing both tools via `setuptools` to ``/home/$user/.local/bin/``

The primer design module as well as the exon analysis and circRNA testing module require a working installation of `R <https://cran.r-project.org/>`_ with `BioConductor <https://www.bioconductor.org/install/>`_. All R packages required are automatically installed during the setup.

.. important:: The setup scripts assumes that the folder for R plugins is writeable (either in the user's home or the system folder).

Python packages
^^^^^^^^^^^^^^^
- For circRNA detection
    * HTSeq>=0.6.1
    * pysam>=0.9.0
    * numpy>=1.8.2
    * pandas>=0.18.1


- For circRNA reconstruction
    * HTSeq>=0.6.1
    * pysam>=0.9.0
    * numpy>=1.8.2
    * pathos>=0.2.1

- For circRNA enrichment
    * pybedtools>=0.7.10
    * statsmodels>=0.8.0

