"""
Provides mathematical functions
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


import copy
import numpy
import math

from pyms.Utils.Error import pymsError
from pyms.Utils.Utils import is_list, is_number


def median(v):
    """
    Returns a median of a list or numpy array

    :param v: Input list or array
    :type v: list or numpy.core.ndarray
    :return: The median of the input list
    :rtype: list

    :author: Vladimir Likic
    """

    if not is_list(v):
        raise TypeError("'v' must be a list or array")

    local_data = copy.deepcopy(v)
    local_data.sort()
    N = len(local_data)

    if (N % 2) == 0:
        # even number of points
        K = N//2 - 1 
        median = (local_data[K] + local_data[K+1])/2.0
    else:
        # odd number of points
        K = (N - 1)//2 - 1
        median = local_data[K+1]

    return median


def vector_by_step(vstart, vstop, vstep):
    """
    generates a list by using start, stop, and step values

    :param vstart: Initial value
    :type vstart: A number
    :param vstop: Max value
    :type vstop: A number
    :param vstep: Step
    :type vstep: A number
   
    :return: A list generated
    :rtype: list

    :author: Vladimir Likic
    """

    if not is_number(vstart) or not is_number(vstop) or not is_number(vstep):
        raise TypeError("parameters 'start', 'stop', and 'step' must be numbers")

    v = []

    p = vstart 
    while p < vstop:
        v.append(p)
        p = p + vstep

    return v


def MAD(v):
    """
    median absolute deviation

    :param v: A list or array
    :type v: list, tuple, or numpy.core.ndarray

    :return: median absolute deviation
    :rtype: float

    :author: Vladimir Likic
    """

    if not is_list(v):
        raise TypeError("'v' must be a list or array")

    m = median(v)
    m_list = []

    for xi in v:
        d = math.fabs(xi - m)
        m_list.append(d)

    mad = median(m_list)/0.6745

    return mad


def amin(v):

    """
    Finds the minimum element in a list or array

    :param v: A list or array
    :type v: list, tuple, or numpy.core.ndarray

    :return: Tuple (maxi, maxv), where maxv is the minimum
        element in the list and maxi is its index
    :rtype: tuple

    :author: Vladimir Likic
    """

    if not is_list(v):
        raise TypeError("'v' must be a list or array")

    minv = max(v) # built-in max() function
    mini = None

    for ii in range(len(v)):
        if v[ii] < minv:
            minv = v[ii]
            mini = ii

    if mini == None:
        raise pymsError("finding maximum failed")

    return mini, minv


def mean(v):

    """
    Calculates the mean

    :param v: A list or array
    :type v: list, tuple, or numpy.core.ndarray

    :return: Mean
    :rtype: float

    :author: Vladimir Likic
    """

    if not is_list(v):
        raise TypeError("'v' must be a list or array")
        
    s = 0.0
    for e in v:
        s = s + e 
    s_mean = s/float(len(v))

    return s_mean


def std(v):

    """
    Calculates standard deviation

    :param v: A list or array
    :type v: list, tuple, or numpy.core.ndarray

    :return: Mean
    :rtype: float

    :author: Vladimir Likic
    """

    if not is_list(v):
        raise TypeError("'v' must be a list or array")

    v_mean = mean(v)

    s = 0.0 
    for e in v:
        d = e - v_mean
        s = s + d*d
    s_mean = s/float(len(v)-1)
    v_std = math.sqrt(s_mean)

    return v_std


def rmsd(list1, list2):

    """
    Calculates RMSD for the 2 lists

    :param list1: First data set
    :type list1: list, tuple, or numpy.core.ndarray
    :param list2: Second data set
    :type list2: list, tuple, or numpy.core.ndarray
    :return: RMSD value
    :rtype: float

    :author: Qiao Wang
    :author: Andrew Isaac
    :author: Vladimir Likic
    """

    if not is_list(list1):
        raise TypeError("'list1' must be a list or array")

    if not is_list(list2):
        raise TypeError("'list2' must be a list or array")

    sum = 0.0
    for i in range(len(list1)):
        sum = sum + (list1[i] - list2[i]) ** 2
    rmsd = math.sqrt(sum / len(list1))
    return rmsd


# added by DK. courtesy of
# http://stackoverflow.com/questions/22354094/pythonic-way-of-detecting-outliers-in-one-dimensional-observation-data
def mad_based_outlier(data, thresh=3.5):
    """
    
    :param data:
    :param thresh:
    :return:
    """
    
    data = numpy.array(data)
    if len(data.shape) == 1:
        data = data[:, None]
    median = numpy.nanmedian(data)
    diff = numpy.nansum((data - median) ** 2, dtype=float, axis=-1)
    diff = numpy.sqrt(diff)
    med_abs_deviation = numpy.nanmedian(diff)
    
    modified_z_score = 0.6745 * diff / med_abs_deviation
    
    return modified_z_score > thresh


# added by DK. courtesy of
# http://stackoverflow.com/questions/22354094/pythonic-way-of-detecting-outliers-in-one-dimensional-observation-data
def percentile_based_outlier(data, threshold=95):
    """
    
    :param data:
    :param threshold:
    :return:
    """
    
    data = numpy.array(data)
    diff = (100 - threshold) / 2.0
    # nanpercentile only works in numpy 1.9 and up
    # minval, maxval = numpy.nanpercentile(data, [diff, 100 - diff])
    data = numpy.array(data)
    minval, maxval = numpy.percentile(numpy.compress(numpy.isnan(data) == False, data), (diff, 100 - diff))
    return (data < minval) | (data > maxval)


# added by DK. courtesy of
# http://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list
def median_outliers(data, m=2.5):
    """
    
    :param data:
    :param m:
    :return:
    """
    
    data = numpy.array(data)
    d = numpy.abs(data - numpy.nanmedian(data))
    mdev = numpy.nanmedian(d)
    s = d / mdev if mdev else 0.
    return (s > m)