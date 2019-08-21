"""
General I/O functions
"""

# Patched version by Dominic Davis-Foster, 2019

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


import types, os, string, sys

import pickle

from pyms.Utils.Utils import is_number, is_str, is_list
from pyms.Utils.Error import pymsError

def dump_object(object, file_name):

    """
    Dumps an object to a file through cPickle.dump()

    :param object: Object to be dumpted
    :type object: An instance of an arbitrary class

    :param file_name: Name of the file for the object dump
    :type file_name: StringType

    :author: Vladimir Likic
    :author: Dominic Davis-Foster
    """

    if not os.path.exists(os.path.dirname(file_name)):
            os.mkdir(os.path.dirname(file_name))
    fp = open_for_writing(file_name, 'wb')
    pickle.dump(object, fp)
    close_for_writing(fp)

def load_object(file_name):

    """
    Loads an object previously dumped with dump_object()

    :param file_name: Name of the object dump file
    :type file_name: StringType

    :return: Object contained in the file 'file_name'
    :rtype: An instance of an arbitrary class

    :author: Vladimir Likic
    """

    fp = open_for_reading(file_name)
    object = pickle.load(fp)
    close_for_reading(fp)

    return object

def open_for_reading(file_name):

    """
    Opens file for reading, returns file pointer

    :param file_name: Name of the file to be opened for reading
    :type file_name: StringType

    :return: Pointer to the opened file
    :rtype: FileType

    :author: Vladimir Likic
    """

    if not is_str(file_name):
        raise TypeError("'file_name' must be a string")
    try:
        fp = open(file_name)
    except IOError:
        raise FileNotFoundError("'%s' does not exist" % (file_name))

    return fp

def open_for_writing(file_name, mode='w'):

    """
    Opens file for writing, returns file pointer

    :param file_name: Name of the file to be opened for writing
    :type file_name: StringType

    :return: Pointer to the opened file
    :rtype: FileType

    :author: Vladimir Likic
    :author: Dominic Davis-Foster
    """

    if not is_str(file_name):
        raise TypeError("'file_name' must be a string")
    #try:
    if not os.path.exists(os.path.dirname(file_name)):
        os.mkdir(os.path.dirname(file_name))
    fp = open(file_name, mode)
#except IOError:
    #    error("Cannot open '%s' for writing" % (file_name))

    return fp

def close_for_reading(fp):

    """
    Closes file pointer open for reading

    :param fp: A file pointer, previously opened for reading
    :type fp: FileType

    :return: none
    :rtype: NoneType

    :author: Vladimir Likic
    """

    fp.close()

def close_for_writing(fp):

    """
    Closes file pointer open for writing

    :param fp: A file pointer, previously opened for writing
    :type fp: FileType

    :return: none
    :rtype: NoneType

    :author: Vladimir Likic
    """

    fp.close()

def file_lines(file_name, filter=False):

    """
    Returns lines from a file, as a list

    :param file_name: Name of a file
    :type: StringType
    :param filter: If True, lines are pre-processes. Newline character
        if removed, leading and taling whitespaces are removed, and lines
        starting with '#' are discarded
    :type: BooleanType

    :return: A list of lines
    :rtype: ListType

    :author: Vladimir Likic
    """

    if not is_str(file_name):
        raise TypeError("'file_name' must be a string")

    fp = open_for_reading(file_name)
    lines = fp.readlines()
    close_for_reading(fp)

    if filter:
        # strip leading and talining whitespaces
        lines_filtered = []
        for line in lines:
            line = line.strip()
            lines_filtered.append(line)

        # discard comments
        lines_to_discard = []
        for line in lines_filtered:
            # remove empty lines and comments
            if len(line) == 0 or line[0] == "#":
                lines_to_discard.append(line)
        for line in lines_to_discard:
            lines_filtered.remove(line)
        lines = lines_filtered

    return lines

def save_data(file_name, data, format_str="%.6f", prepend="", sep=" ",
    compressed=False):

    """
    Saves a list of numbers or a list of lists of numbers
        to a file with specific formatting

    :param file_name: Name of a file
    :type: StringType
    :param data: A list of numbers, or a list of lists
    :type: ListType
    :param format_str: A format string for individual entries
    :type: StringType
    :param prepend: A string, printed before each row
    :type: StringType
    :param sep: A string, printed after each number
    :type: StringType
    :param compressed: A boolean. If True, the output will be gzipped
    :type: BooleanType

    :return: none
    :rtype: NoneType

    :author: Vladimir Likic
    :author: Dominic Davis-Foster
    """

    if not os.path.exists(os.path.dirname(file_name)):
            os.mkdir(os.path.dirname(file_name))
    
    if not is_str(file_name):
        raise TypeError("'file_name' must be a string")

    if not is_list(data):
        raise TypeError("'data' must be a list")

    if not is_str(prepend):
        raise TypeError("'prepend' must be a string")

    if not is_str(sep):
        raise TypeError("'sep' must be a string")

    fp = open_for_writing(file_name)

    # decide whether data is a vector or matrix
    if is_number(data[0]):
        for item in data:
            if not is_number(item):
                raise TypeError("not all elements of the list are numbers")
        data_is_matrix = 0
    else:
        for item in data:
            if not is_list(item):
                raise TypeError("not all elements of the list are lists")
        data_is_matrix = 1

    if data_is_matrix:
        for ii in range(len(data)):
            fp.write(prepend)
            for jj in range(len(data[ii])):
                if is_number(data[ii][jj]):
                    fp.write(format_str % (data[ii][jj]))
                    if (jj<(len(data[ii])-1)): fp.write(sep)
                else:
                    raise TypeError("'datum' must be a number")
            fp.write("\n")
    else:
        for ii in range(len(data)):
            fp.write(prepend)
            fp.write(format_str % (data[ii]))
            fp.write("\n")

    close_for_writing(fp)

    if compressed:
        status = os.system('gzip %s' % (file_name))
        if status != 0:
            pymsError("gzip compress failed")

