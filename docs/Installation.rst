Installation
********************************************************


**circtools** is written in Python2 (>= 2.7, ``detect`` and ``reconstruct`` module) and Python3 (>=3.4, all other modules). The tool has a number of external dependencies, mostly standard bioinformatics tools and packages. The installation will, by default, try to install all required dependencies.

Installation is performed via `python3 setup.py install`. No sudo access is required if the installation is suffixed with ``--user`` which will install the package in a user-writeable folder. In this case, the binaries should be installed to ``/home/$USER/.local/bin/`` (for Debian-based systems).

``circtools`` was developed and tested on Debian Jessie 8 64 Bit. macOS support is currently (08/2018) being tested and is already available in the ``mac-dev`` branch of the github repository (however, the macOS functionality cannot be fully guaranteed yet). 

Installation from PyPi (preferred)
-----------------------------------

The default installation will install everything needed to run circtools *except R, STAR, or Stringtie* (see below). If you like you may install circtools locally (first call) or globally (second call, SU required).

.. code-block:: bash

    pip3 install circtools --user # does not require root access, installation to local user directory

.. code-block:: bash

    pip3 install circtools # will require root access and globally install circtools


Installation from GitHub
--------------------------

The GitHub installation will install the most recent version directly from the source repository. Use this method if you want the latest fixes and features.

.. code-block:: bash

    git clone https://github.com/dieterich-lab/circtools.git
    cd circtools
    pip3 install . --verbose --user

Updating circtools
--------------------------

You may want to update the circtools package if new versions are published. Like for the initial installation there are two ways to update circtools:

.. code-block:: bash

    pip3 install circtools --user --upgrade

.. code-block:: bash

    cd /path/to/circtools/repo/
    git pull
    pip3 install . install --verbose --user --upgrade

Required dependencies
---------------------

External tools
^^^^^^^^^^^^^^^

* `bedtools [>= 2.27.1] <http://bedtools.readthedocs.io/en/latest/content/installation.html>`_ required by the enrichment module

* `R [>= 3.3] <https://www.digitalocean.com/community/tutorials/how-to-install-r-on-ubuntu-16-04-2>`_ required by visualisation scripts and the primer design module

* `STAR [>= 2.6.0] <https://github.com/alexdobin/STAR>`_ required by the ``detect`` and ``reconstruct`` module to map RNA-seq reads against a reference genome and detect back splice junctions

* `Stringtie [>= 1.3.3b, optional] <https://github.com/gpertea/stringtie>`_ required by the ``exon`` module to carry out exon level analyses.

The installation procedure will automatically install two additional Python-based dependencies: `DCC <https://github.com/dieterich-lab/DCC>`_ and `FUCHS <https://github.com/dieterich-lab/FUCHS>`_ by temporarily cloning the repositories and installing both tools via `setuptools` to ``/home/$USER/.local/bin/``. Both tools **require Python 2** in order to run.

The primer design module as well as the exon analysis and circRNA testing module require a working installation of `R <https://cran.r-project.org/>`_ with `BioConductor <https://www.bioconductor.org/install/>`_. All R packages required are automatically installed during the setup.

.. important:: The setup scripts assumes that the folder for R plugins is writeable (either in the user's home or the system folder).

Python packages
^^^^^^^^^^^^^^^
- For circRNA detection
    * HTSeq>=0.11.0
    * pysam>=0.13.0
    * numpy>=1.8.2
    * pandas>=0.18.1

- For circRNA reconstruction
    * HTSeq>=0.11.0
    * pysam>=0.13.0
    * numpy>=1.8.2
    * pathos>=0.2.1

- For circRNA enrichment
    * pybedtools>=0.7.10
    * statsmodels>=0.8.0

- For circRNA primer design
    * BioPython>=1.71


Detailed installation
----------------------

Getting the source code
^^^^^^^^^^^^^^^^^^^^^^^

**Step 1**: Clone source code from GitHub:

.. code-block:: bash

    git clone https://github.com/dieterich-lab/circtools.git

Installation
^^^^^^^^^^^^

**Step 2**: Install circtools using the provided installation script. The ``--user`` flag installs circtools in your home folder, thus making sure you do not require any administrative rights during the installation:

.. code-block:: bash

    cd circtools
    pip3 install . install --verbose --user

R environment
^^^^^^^^^^^^^^

**Step 3**: Setting up R environment. In order for the automatic installation of R packages to work we need to set the package directory to a user-writeable path. The setup automatically sets that path to ``/home/$USER/.R/``.


Dependencies
^^^^^^^^^^^^

**Step 4**: The setup script is designed to make sure that the environment is setup correctly to run circtools. The circtools setup will automatically install `CircTest <https://github.com/dieterich-lab/CircTest>`_, `primex <https://github.com/dieterich-lab/primex>`_, `DCC <https://github.com/dieterich-lab/DCC>`_ and `FUCHS <https://github.com/dieterich-lab/FUCHS>`_.

Finishing up
^^^^^^^^^^^^

**Step 5**: Adding installation folder to ``$PATH``. In order for circtools to find all executables, the setup will add the folder ``/home/$USER/.local/bin/`` automatically to your ``.bashrc`` file

This closes the circtools installation. To verify that circtools has been correctly installed, try to call circtools for the first time:

.. code-block:: bash

    $> circtools --help
    usage: circtools [-V] <command> [<args>]
