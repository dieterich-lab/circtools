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


class PrimerDesign(circ_module.circ_template.CircTemplate):
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

        # needed for Rscript decoupling
        import subprocess

        # param1 = self.cli_params.output_directory

        try:
            # we assume Rscript is in the user's $PATH variable

            # But..let's check that first, you never know

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

            exit()
            # ------------------------------------ need to call the correct R script here -----------------------

            # need to define path top R wrapper
            path2script = 'path/to your script/max.R'

            # Variable number of args in a list
            args = ['11', '3', '9', '42']

            # Build subprocess command
            cmd = [self.command, path2script] + args
            #
            # check_output will run the command and store to result

            # ------------------------------------ run script and check output -----------------------

            module_output = subprocess.check_output(cmd, universal_newlines=True)

        except subprocess.CalledProcessError:
            self.log_entry("`Rscript` (a part of the R installation) could not be found in your $PATH.")


