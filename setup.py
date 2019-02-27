#! /usr/bin/env python3

"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from os import path
from setuptools.command.develop import develop
from setuptools.command.install import install
from setuptools import setup, find_packages

import sys


def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        print(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            print("Please respond with 'yes' or 'no' "
                  "(or 'y' or 'n').\n")


class PostDevelopCommand(develop):
    """Post-installation for development mode."""

    def run(self):
        develop.run(self)


class PostInstallCommand(install):
    """Post-installation for installation mode."""

    def run(self):

        import subprocess

        import os
        # read the docs does not need any further dependencies, just the doc code
        if not os.environ.get('READTHEDOCS') == 'True':

            subprocess.check_call(["bash", "scripts/install_create_r_environ.sh"])
            subprocess.check_call(["bash", "scripts/install_external.sh"])
            subprocess.check_call(["bash", "scripts/install_add_to_bashrc.sh"])

        install.run(self)
        # place for post install commands


here = path.abspath(path.dirname(__file__))

if sys.version_info.major < 3:
    sys.exit('\n'
             'Sorry, Python 2 is not supported\n'
             'Did you run pip install circtools?\n'
             'Try \'pip3 install circtools\'')

elif sys.version_info.minor < 4:
    sys.exit('\nSorry, Python < 3.4 is not supported\n')

# Get the long description from the relevant file
with open(path.join(here, 'README.rst')) as f:
    long_description = f.read()

setup(
    name='circtools',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version='1.1.0.8',

    description='circtools - a circular RNA toolbox',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/dieterich-lab/circtools',

    # Author details

    author='Tobias Jakobi',
    author_email='Tobias.Jakobi@med.Uni-Heidelberg.DE',

    # Choose your license
    license='GNU General Public License (GPL)',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are

        # Development Status :: 1 - Planning
        # Development Status :: 2 - Pre-Alpha
        # Development Status :: 3 - Alpha
        # Development Status :: 4 - Beta
        # 'Development Status :: 5 - Production/Stable',
        # Development Status :: 6 - Mature
        # Development Status :: 7 - Inactive

        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        # 'Programming Language :: Python :: 2',
        # 'Programming Language :: Python :: 2.6',
        # 'Programming Language :: Python :: 2.7',
        # 'Programming Language :: Python :: 3',
        # 'Programming Language :: Python :: 3.2',
        # 'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',

    ],

    # What does your project relate to?
    keywords='circular RNA bioinformatics',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(),

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=[
        'HTSeq>=0.11.0',
        'pysam >= 0.13',
        'numpy>=1.14.5',
        'pybedtools>=0.7.10',
        'biopython >= 1.71',
        'scipy>=0.19.0'
    ],

    python_requires='>=3.4',

    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    # extras_require={
    #     'dev': ['check-manifest'],
    #     'test': ['coverage'],
    # },

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_data={
        'circtools': ['scripts/circtools'],
    },

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    # data_files=[('my_data', ['data/data_file'])],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.

    zip_safe=False,

    # entry_points={
    #     'console_scripts': [
    #         'circtools=circtools:main'
    #     ],
    # },

    # this will be our main "executable"
    scripts=[
        'scripts/circtools',
        'scripts/circtools_detect_write_skip_tracks.pl',
        'scripts/circtools_enrich_visualization.R',
        'scripts/circtools_quickcheck_wrapper.R',
        'scripts/circtools_reconstruct_visualization.R',
        'scripts/circtools_primex_wrapper.R',
        'scripts/circtools_primex_formatter.R',
        'scripts/circtools_exon_wrapper.R',
        'scripts/circtools_circtest_wrapper.R',
        'scripts/circtools_generate_intron_gtf.sh',
        'scripts/circtools_generate_flanking_introns.py',
        'scripts/create_igv_script_from_gene_names.py',
        'scripts/create_igv_script_from_position_list.py',
        'scripts/create_igv_script.py',
        'scripts/detect_new_exons_from_fuchs_data.py',
        'scripts/get_introns_from_ensembl.pl'
    ],

    # adding support for post install scripts from
    # https://stackoverflow.com/questions/20288711/post-install-script-with-python-setuptools

    cmdclass={
        'develop': PostDevelopCommand,
        'install': PostInstallCommand,
    },

    project_urls={  # Optional
        'Bug Reports': 'https://github.com/dieterich-lab/circtools/issues',
        'Dieterich Lab': 'https://dieterichlab.org',
        'Source': 'https://github.com/dieterich-lab/circtools',
        'Documentation': 'http://docs.circ.tools'
},
)
