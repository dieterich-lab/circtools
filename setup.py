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
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")


class PostDevelopCommand(develop):
    """Post-installation for development mode."""

    def run(self):

        develop.run(self)


class PostInstallCommand(install):
    """Post-installation for installation mode."""

    def run(self):

        import subprocess
        from time import sleep

        print("")

        # step 1: install DCC, FUCHS and primer design R module
        # print("We need to install two other programs of the Dieterich lab circRNA suit: DCC and FUCHS.")
        # print("We'll install both of them for you from GitHub.")
        # print("")
        # print("If you want to cancel the installation press CTRL-C now.")
        # sleep(10)
        # import subprocess
        # subprocess.check_call(["sh", "scripts/external_install.sh"])

        # step 2: update $PATH
        print("In order for circtools to be globally callable, we would add the installation folder to the $PATH")
        print("variable. Would you like us to do that?")

        answer = query_yes_no("Update $PATH in .bashrc?")
        if answer:
            print("Running update script...")
            subprocess.check_call(["sh", "scripts/add_to_bashrc.sh"])
        else:
            print("Okay. Please update the $PATH variable yourself, otherwise you may not be able to run circtools.")
        sleep(5)

        print("")

        # step 3: create .Renviron file
        print("In order for the circtools primer design module to run, we need to install some R modules.")
        print("Please make sure R >= 3.4.0 is installed.")
        print("Should we update the R package location in order to install package as user?")

        answer = query_yes_no("Update R_LIB in .Renviron")
        if answer:
            print("Running update script...")
            subprocess.check_call(["sh", "scripts/create_r_environ.sh"])
        else:
            print("Okay. If the R library path is not writeable the installation will most probably fail.")
        sleep(5)

        print("")

        # step 4: install OligoArrayAux
        print("The circtools primer module also requires the software \"OligoArrayAux\".")
        print("Should we install the software for you?")
        print("If you choose yes OligoArrayAux will be installed into $HOME/.local/.")

        answer = query_yes_no("Install OligoArrayAux?")
        if answer:
            print("Running installer script...")
            subprocess.check_call(["sh", "scripts/install_oligo_array_aux.sh"])
        else:
            print("Okay. Please install OligoArrayAux yourself. Otherwise the primer module will not be functional.")
        sleep(5)

        install.run(self)
        # place for post install commands


here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
# with open(path.join(here, 'DESCRIPTION.rst'), encoding='utf-8') as f:
#     long_description = f.read()

setup(
    name='circtools',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version='0.0.1',

    description='circtools - a circular RNA toolbox',
    # long_description=long_description,

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
        # Development Status :: 5 - Production/Stable
        # Development Status :: 6 - Mature
        # Development Status :: 7 - Inactive

        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Researchers',
        'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish (should match "license" above)
        'License :: GNU General Public License (GPL) 3.0',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        # 'Programming Language :: Python :: 2',
        # 'Programming Language :: Python :: 2.6',
        # 'Programming Language :: Python :: 2.7',
        # 'Programming Language :: Python :: 3',
        # 'Programming Language :: Python :: 3.2',
        # 'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
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
        'pybedtools>=0.7.10',
        'statsmodels>=0.8.0'
    ],

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
        'scripts/circtools_primer',
    ],

    # adding support for post install scripts from
    # https://stackoverflow.com/questions/20288711/post-install-script-with-python-setuptools

    cmdclass={
        'develop': PostDevelopCommand,
        'install': PostInstallCommand,
    },
)
