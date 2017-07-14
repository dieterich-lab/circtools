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

import time
import logging
import os
import sys
from abc import ABCMeta, abstractmethod


class CircTemplate(object):
    """A template class for new sub modules.

    Attributes:
        cli_params: Command lines params supplied by the user for the current module
        program_name: Name of this module
        version: Version of the module
    """
    __metaclass__ = ABCMeta

    def __init__(self, argparse_arguments, program_name, version):

        # get the user supplied options
        self.cli_params = argparse_arguments
        self.program_name = program_name
        self.version = version

        # set time format
        time_format = time.strftime("%Y_%m_%d__%H_%M")

    @staticmethod
    def log_entry(string):
        """Logs to log file and prints on screen
        """
        message = string
        logging.info(message)
        print(message)

    @staticmethod
    def check_input_files(input_file_list):
        """Checks supplied list of files for existence.
        Will halt the program if file not accessible
        """
        for file in input_file_list:
            # check if exists
            if not os.path.isfile(file):
                message = ("File " + str(file) + " cannot be found, exiting.")
                logging.info(message)
                sys.exit(message)

    @staticmethod
    def check_int_arguments(input_list):
        """Checks supplied list of files for existence.
        Will halt the program if file not accessible
        """
        for number in input_list:
            # check if exists
            try:
                int(number)
            except ValueError:
                message = ("Error: column %s is no valid column index." % str(number))
                logging.info(message)
                sys.exit(message)

    @staticmethod
    def check_float_arguments(input_list):
        """Checks supplied list of files for existence.
        Will halt the program if file not accessible
        """
        for number in input_list:
            # check if exists
            try:
                float(number)
            except ValueError:
                message = ("Error: column %s is no valid column index." % str(number))
                logging.info(message)
                sys.exit(message)

    @abstractmethod
    def module_name(self):
        """"Return a string representing the type of vehicle this is."""
        pass

    @abstractmethod
    def run_module(self):
        pass
