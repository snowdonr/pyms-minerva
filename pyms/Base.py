"""
Base for PyMassSpec classes.
"""

################################################################################
#                                                                              #
#    PyMassSpec software for processing of mass-spectrometry data              #
#    Copyright (C) 2005-2012 Vladimir Likic                                    #
#    Copyright (C) 2019-2020 Dominic Davis-Foster                              #
#                                                                              #
#    This program is free software; you can redistribute it and/or modify      #
#    it under the terms of the GNU General Public License version 2 as         #
#    published by the Free Software Foundation.                                #
#                                                                              #
#    This program is distributed in the hope that it will be useful,           #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
#    GNU General Public License for more details.                              #
#                                                                              #
#    You should have received a copy of the GNU General Public License         #
#    along with this program; if not, write to the Free Software               #
#    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 #
#                                                                              #
################################################################################

# stdlib
import typing
import pathlib
import os

# this package
from pyms.Utils.IO import prepare_filepath
from pyms.Utils.Utils import _pickle_dump_path, is_path

__all__ = ["pymsBaseClass"]
PathLike = typing.Union[str, pathlib.Path, os.PathLike]


class pymsBaseClass:
    """
    Base class.
    """

    def dump(self, file_name: PathLike, protocol: int = 3):
        """
        Dumps an object to a file through :func:`pickle.dump()`.

        :param file_name: Filename to save the dump as.
        :param protocol: The pickle protocol to use.

        :authors: Vladimir Likic, Dominic Davis-Foster (pathlib and pickle protocol support)
        """  # noqa: D402  # TODO: False positive

        if not is_path(file_name):
            raise TypeError("'file_name' must be a string or a PathLike object")

        file_name = prepare_filepath(file_name)
        _pickle_dump_path(file_name, self, protocol)
