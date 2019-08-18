"""
Functions related to storing and loading a list of Peak objects
"""

#############################################################################
#                                                                           #
#    PyMS software for processing of metabolomic mass-spectrometry data     #
#    Copyright (C) 2005-2012 Vladimir Likic                                 #
#    Copyright (C) 2019 Dominic Davis-Foster                                #
#                                                                           #
#    This program is free software; you can redistribute it and/or modify   #
#    it under the terms of the GNU General Public License version 2 as      #
#    published by the Free Software Foundation.                             #
#                                                                           #
#    This program is distributed in the hope that it will be useful,        #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of         #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
#    GNU General Public License for more details.                           #
#                                                                           #
#    You should have received a copy of the GNU General Public License      #
#    along with this program; if not, write to the Free Software            #
#    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.              #
#                                                                           #
#############################################################################


import pickle
import pathlib

from pyms.Peak.Class import Peak
from pyms.Utils.Utils import is_list


def store_peaks(peak_list, file_name):
    """
    :summary:Store the list of peak objects

    :param peak_list: A list of peak objects
    :type peak_list: pyms.Peaks.Class.Peak
    :param file_name: File name to store peak list
    :type file_name: str or pathlib.Path

    :author: Andrew Isaac
    :author: Dominic Davis-Foster (pathlib support)
    """

    if not isinstance(file_name, (str, pathlib.Path)):
        raise TypeError("'file_name' must be a string or a pathlib.Path object")

    if not isinstance(file_name, pathlib.Path):
        file_name = pathlib.Path(file_name)

    fp = file_name.open('wb')
    pickle.dump(peak_list, fp, 1)
    fp.close()


def load_peaks(file_name):

    """
    :summary: Loads the peak_list stored with 'store_peaks'

    :param file_name: File name of peak list
    :type file_name: str or pathlib.Path


    :return: The list of Peak objects
    :rtype: list

    :author: Andrew Isaac
    :author: Dominic Davis-Foster (pathlib support)
    """

    if not isinstance(file_name, (str, pathlib.Path)):
        raise TypeError("'file_name' must be a string or a pathlib.Path object")

    if not isinstance(file_name, pathlib.Path):
        file_name = pathlib.Path(file_name)

    fp = file_name.open('rb')
    peak_list = pickle.load(fp)
    fp.close()

    if not is_list(peak_list):
        raise TypeError("'file_name' is not a List")
    if not len(peak_list) > 0 and not isinstance(peak_list[0], Peak):
        raise TypeError("'peak_list' must be a list of Peak objects")

    return peak_list
