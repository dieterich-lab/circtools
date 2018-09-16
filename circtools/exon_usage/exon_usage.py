#! /usr/bin/env python3

# Copyright (C) 2017 Tobias Jakobi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either self.version 3 of the License, or
# (at your option) any later self.version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import circ_module.circ_template


class ExonUsage(circ_module.circ_template.CircTemplate):
    def __init__(self, argparse_arguments, program_name, version):

        # get the user supplied options
        self.cli_params = argparse_arguments
        self.program_name = program_name
        self.version = version
        self.command = 'Rscript'

    def module_name(self):
        """"Return a string representing the name of the module."""
        return self.program_name

    def run_module(self):

        # we first need to make sure all parameters are sane, the R script does not validate any parameters

        # check output directory

        import os

        # let's first check if the output directory exists
        if not (os.access(self.cli_params.output_directory, os.W_OK)):
            self.log_entry("Output directory %s not writable." % self.cli_params.output_directory)
            # exit with -1 error if we can't use it
            exit(-1)

        # let's first check if the temporary directory exists
        if not os.path.exists(self.cli_params.output_directory):
            os.makedirs(self.cli_params.output_directory)

        # check DCC directory
        if not (os.path.exists(self.cli_params.DCC_dir)):
            self.log_entry("DCC/detect data directory %s does not exist (or is not accessible)."
                           % self.cli_params.DCC_dir)
            # exit with -1 error if we can't use it
            exit(-1)

        # check DCC files (only existence, not the content)
        self.check_input_files([
                                self.cli_params.DCC_dir+"CircRNACount",
                                self.cli_params.DCC_dir+"LinearCount",
                                self.cli_params.DCC_dir+"CircCoordinates"
                                ])

        # check Ballgown directory
        if self.cli_params.ballgown_data and not (os.path.exists(self.cli_params.ballgown_data)):
            self.log_entry("Ballgown data directory %s does not exist (or is not accessible)."
                           % self.cli_params.ballgown_data)
            # exit with -1 error if we can't use it
            exit(-1)

        # check circtest and GTF files (only existence, not the content)
        self.check_input_files([self.cli_params.gtf_file,  self.cli_params.circtest_file])

        # check sample names
        if len(self.cli_params.condition_list.split(",")) < 2:
            self.log_entry("Error: Length of parameter list specified via -c is < 2.")

        # check columns
        if len(self.cli_params.condition_columns.split(",")) < 2:
            self.log_entry("Error: Length of parameter list specified via -l is < 2.")

        # check grouping
        if len(self.cli_params.grouping.split(",")) < 2:
            self.log_entry("Error: Length of parameter list specified via -g is < 2.")

        for column in self.cli_params.condition_columns.split(","):
            try:
                int(column)
            except ValueError:
                self.log_entry("Error: column %s is no valid column index." % str(column))
                exit(-1)

        for column in self.cli_params.grouping.split(","):
            try:
                int(column)
            except ValueError:
                self.log_entry("Error: group %s is no valid group index." % str(column))
                exit(-1)

        # check numeric arguments

        # self.check_int_arguments([self.cli_params.max_plots])

        # needed for Rscript decoupling
        import subprocess

        # import re module
        import re

        r_location = subprocess.check_output(['which', self.command], universal_newlines=True,
                                             stderr=subprocess.STDOUT).split('\n')[0]

        r_version = subprocess.check_output([self.command, '--version'], universal_newlines=True,
                                            stderr=subprocess.STDOUT)
        # okay, Rscript is really there, we put together the command line now:

        m = re.search('(\d+\.\d+\.\d+)', r_version)
        r_version = m.group(0)

        self.log_entry("Using R version %s [%s]" % (r_version, r_location))

        # ------------------------------------ need to call the correct R script here -----------------------

        # need to define path top R wrapper
        exon_script = 'circtools_exon_wrapper.R'

        # Variable number of args in a list
        args = [
                self.cli_params.DCC_dir,
                self.cli_params.replicates,
                self.cli_params.condition_list,
                self.cli_params.condition_columns,
                self.cli_params.grouping,
                self.cli_params.output_directory+"/"+self.cli_params.output_prefix,
                self.cli_params.ballgown_data,
                self.cli_params.gtf_file,
                self.cli_params.circtest_file,
                self.cli_params.has_header,
                self.cli_params.species
        ]

        # ------------------------------------ run script and check output -----------------------

        import os
        os.system(exon_script + " " + ' '.join(str(e) for e in args))
