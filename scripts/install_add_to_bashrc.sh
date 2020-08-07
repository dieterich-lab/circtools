#!/bin/bash

# Copyright (C) 2017 Tobias Jakobi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# we append the local user python path to the $PATH variable
# this make circtools callable for the user

# are we running in an virtual environment?

CIRCTOOLS=`which circtools`

if [ ! "$CIRCTOOLS" ]; then

    if [ $VIRTUAL_ENV ]; then
        echo "export PATH=\$PATH:$VIRTUAL_ENV/bin" >> ~/.bashrc
    else
        echo "export PATH=\$PATH:~/.local/bin" >> ~/.bashrc
    fi
fi