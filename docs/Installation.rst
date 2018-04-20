Installation
=============


**circtools** is written in Python3 (>=3.4). The tool has a number of external dependencies, mostly standard bioinformatics tools and packages. The installation will, by default, try to install all required dependencies.

Installation is performed via `python3 setup.py install`. No sudo access is required if the installation is suffixed with ``--user`` which will install the package in a user-writeable folder. In this case, the binaries should be installed to ``/home/$USER/.local/bin/`` (for Debian-based systems).

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


* `bedtools [>= 2.27.1] <http://bedtools.readthedocs.io/en/latest/content/installation.html>`_ required by the enrichment module

.. important:: Currently at least bedtools v2.26.0-125-g52db654-dirty is required to run the enrichment module. Earlier versions have issues with the ``shuffle`` command.

* `R [>= 3.3] <https://www.digitalocean.com/community/tutorials/how-to-install-r-on-ubuntu-16-04-2>`_ required by visualisation scripts and the primer design module


The installation procedure will automatically install two additional Python-based dependencies: `DCC <https://github.com/dieterich-lab/DCC>`_ and `FUCHS <https://github.com/dieterich-lab/FUCHS>`_ by temporarily cloning the repositories and installing both tools via `setuptools` to ``/home/$USER/.local/bin/``

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

- For circRNA detection
    * BioPython>=1.71


Detailed installation
----------------------

Getting the source code
^^^^^^^^^^^^^^^^^^^^^^^

**Step 1**: Clone source code from GitHub::

    git clone git@github.com:dieterich-lab/circtools.git

**Step 2**: Install circtools using the provided installation script. The ``--user`` flag installs circtools in your home folder, thus making sure you do not require any administrative rights during the installation::

    cd circtools
    python3 setup.py install --verbose --user


**Step 3**: Setting up R environment. In order for the automatic installation of R packages to work we need to set the package directory to a user-writeable path. The setup can automatically set that path to /home/$USER/.R/::

    Should we update the R package location in order to install package as user?
    Update R_LIB in .Renviron [Y/n]

**Step 4**: The setup script is designed to guide you through the installation process and makes sure your enviroonment is setup correctly to run circtools. You will have to answer a few questions throughout this process::

    We need to install two other programs of the Dieterich Lab circRNA suit, DCC and FUCHS, as well as R package dependencies for other modules of circtools
    We'll install everything for you from GitHub and CRAN for you.
    
    In order for the circtools primer design module to run, we need to install some R modules.
    Please make sure R >= 3.3 is installed and your R library path is writeable .
    
    Do you want to continue the automatic dependency installation?
    -> "n" will only install the circtools base package
    -> CTRL-C will abort the installation
     [Y/n]

Answer with "y" to automatically install `CircTest <https://github.com/dieterich-lab/CircTest>`_, `primex <https://github.com/dieterich-lab/primex>`_, `DCC <https://github.com/dieterich-lab/DCC>`_ and `FUCHS <https://github.com/dieterich-lab/FUCHS>`_. 

**Step 5**: Adding installation folder to $PATH. In order for circtools to find all exectuables, the setup will give you the possibility to add the folder ``/home/$USER/.local/bin/`` automatically to your ``.bashrc`` file::

    In order for circtools to be globally callable, we would add the installation folder to the $PATH variable. Would you like us to do that?
    Update $PATH in .bashrc? [Y/n]

This closes the circtools installation. To verify that circtools has been correctly installed, try to call circtools for the first time::

    $> circtools --help
    usage: circtools [-V] <command> [<args>]


