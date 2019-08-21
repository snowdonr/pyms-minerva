"""
Classes to model a GC-MS Ion Chromatogram
"""

#############################################################################
#                                                                           #
#    PyMassSpec software for processing of metabolomic mass-spectrometry data     #
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


import numpy
import math
import copy

import deprecation
import warnings
from pyms import __version__

from pyms.Utils.Error import pymsError
from pyms.Utils.Utils import is_list
from pyms.Utils.IO import open_for_writing, close_for_writing


class IonChromatogram(object):
	"""
	Models an ion chromatogram

		An ion chromatogram is a set of intensities as a function of retention
		time. This can can be either m/z channel intensities (for example, ion
		chromatograms at m/z=65), or cumulative intensities over all measured
		m/z. In the latter case the ion chromatogram is total ion chromatogram
		(TIC).

		The nature of an IonChromatogram object can be revealed by inspecting
		the value of the attribute '__mass'. This is set to the m/z value of the
		ion chromatogram, or to None for TIC.

	:author: Lewis Lee
	:author: Vladimir Likic
	:author: Dominic Davis-Foster (type assertions and properties)
	"""
	
	def __init__(self, ia, time_list, mass=None):
		"""
		:param ia: Ion chromatogram intensity values
		:type ia: numpy.array
		:param time_list: A list of ion chromatogram retention times
		:type time_list: list
		:param mass: Mass of ion chromatogram (Null if TIC)
		:type mass: int or float

		:author: Lewis Lee
		:author: Vladimir Likic
		"""
		
		if not isinstance(ia, numpy.ndarray):
			raise TypeError("'ia' must be a numpy array")
		
		if not is_list(time_list) or not isinstance(time_list[0], (int, float)):
			raise TypeError("'time_list' must be a list of numbers")
		
		if len(ia) != len(time_list):
			raise ValueError("Intensity array and time list differ in length")
		
		if mass and not isinstance(mass, (int, float)):
			raise TypeError("'mass' must be a number")
		
		self.__ia = ia
		self.__time_list = time_list
		self.__mass = mass
		self.__time_step = self.__calc_time_step(time_list)
		self.__min_rt = min(time_list)
		self.__max_rt = max(time_list)
	
	def __len__(self):
		
		"""
		Returns the length of the IonChromatogram object

		:return: Length of ion chromatogram
		:rtype: int

		:author: Lewis Lee
		:author: Vladimir Likic
		"""
		
		return self.__ia.size
	
	def __sub__(self, other):
		"""
		Subtracts another IC from the current one

		:param other: Another IC
		:type other: pyms.GCMS.IonChromatogram
		"""
		
		ia_for_sub = other.intensity_array
		
		for i in range(self.__ia.size):
			self.__ia[i] = self.__ia[i] - ia_for_sub[i]
		
		return self
	
	def __eq__(self, other):
		if isinstance(other, self.__class__):
			return self.time_list == other.time_list \
				   and all(numpy.equal(self.intensity_array, other.intensity_array)) \
				   and self.mass == other.mass
			
		return NotImplemented
	
	def __copy__(self):
		return IonChromatogram(ia=numpy.copy(self.__ia), time_list=self.__time_list[:], mass=copy.copy(self.__mass))
	
	def __deepcopy__(self, memodict={}):
		return self.__copy__()
	
	def get_index_at_time(self, time):
		"""
		Returns the nearest index corresponding to the given time

		:param time: Time in seconds
		:type time: float

		:return: Nearest index corresponding to given time
		:rtype: int

		:author: Lewis Lee
		:author: Tim Erwin
		:author: Milica Ng
		:author: Vladimir Likic
		"""
		
		if not isinstance(time, (int, float)):
			raise TypeError("'time' must be a number")
		
		if time < self.__min_rt or time > self.__max_rt:
			raise IndexError("time %.2f is out of bounds (min: %.2f, max: %.2f)" %
							 (time, self.__min_rt, self.__max_rt))
		
		time_list = self.__time_list
		time_diff_min = self.__max_rt
		ix_match = None
		
		for ix in range(len(time_list)):
			
			time_diff = math.fabs(time - time_list[ix])
			
			if time_diff < time_diff_min:
				ix_match = ix
				time_diff_min = time_diff
		
		return ix_match
	
	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'IonChromatogram.intensity_array' instead")
	def get_intensity_array(self):
		"""
		Returns the entire intensity array

		:return: Intensity array
		:rtype: numpy.ndarray

		:author: Lewis Lee
		:author: Vladimir Likic
		"""
		
		return self.intensity_array
		
	def get_intensity_at_index(self, ix):
		"""
		Returns intensity at given index

		:param ix: An index
		:type ix: int

		:return: Intensity value
		:rtype: float

		:author: Lewis Lee
		:author: Vladimir Likic
		"""
		
		if not isinstance(ix, int):
			raise TypeError("'ix' must be an integer")
		
		if ix < 0 or ix > self.__ia.size - 1:
			raise IndexError("index out of bounds")
		
		return self.__ia[ix]
	
	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'IonChromatogram.mass' instead")
	def get_mass(self):
		"""
		Returns the m/z channel of the IC

		:return: m/z channel of the IC
		:rtype: int

		:author: Sean O'Callaghan
		"""
		if self.__mass == None:
			raise pymsError("TIC has no m/z label")
		
		return self.mass
	
	def get_time_at_index(self, ix):
		"""
		Returns time at given index

		:param ix: An index
		:type ix: int

		:return: Time value
		:rtype: float

		:author: Lewis Lee
		:author: Vladimir Likic
		"""
		
		if not isinstance(ix, (int, float)):
			raise TypeError("'ix' must be an integer")
		
		if ix < 0 or ix > len(self.__time_list) - 1:
			raise IndexError("index out of bounds")
		
		return self.__time_list[ix]
	
	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'IonChromatogram.time_list' instead")
	def get_time_list(self):
		"""
		Returns the time list

		:return: Time list
		:rtype: list

		:author: Lewis Lee
		:author: Vladimir Likic
		"""
		
		return self.time_list
	
	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'IonChromatogram.time_step' instead")
	def get_time_step(self):
		"""
		Returns the time step

		:return: Time step
		:rtype: float

		:author: Lewis Lee
		:author: Vladimir Likic
		"""
		
		return self.__time_step
		
	@property
	def intensity_array(self):
		"""
		Returns the entire intensity array

		:return: Intensity array
		:rtype: numpy.ndarray

		:author: Lewis Lee
		:author: Vladimir Likic
		"""
		
		return numpy.copy(self.__ia)
	
	@intensity_array.setter
	def intensity_array(self, ia):
		"""
		Sets the value for the intensity array

		:param ia: An array of new intensity values
		:type ia: numpy.ndarray

		:return: none
		:rtype: NoneType

		:author: Vladimir Likic
		"""
		
		self.__ia = ia
	
	def is_tic(self):
		"""
		Returns True if the ion chromatogram is a total ion
			chromatogram (TIC), or False otherwise

		:return: A boolean value indicating if the ion chromatogram
			is a total ion chromatogram (True) or not (False)
		:rtype: bool

		:author: Lewis Lee
		:author: Vladimir Likic
		"""
		
		return self.__mass is None
	
	@property
	def mass(self):
		"""
		Returns the m/z channel of the IC

		:return: m/z channel of the IC
		:rtype: int or float

		:author: Sean O'Callaghan
		"""
		if self.__mass is None:
			warnings.warn("TIC has no m/z label", Warning)
		
		return self.__mass
	
	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'IonChromatogram.intensity_array' instead")
	def set_intensity_array(self, ia):
		"""
		Sets the value for the intensity array

		:param ia: An array of new intensity values
		:type ia: numpy.ndarray

		:return: none
		:rtype: NoneType

		:author: Vladimir Likic
		"""
		# todo: type check
		self.__ia = ia
		
	@property
	def time_list(self):
		"""
		Returns the time list

		:return: Time list
		:rtype: list

		:author: Lewis Lee
		:author: Vladimir Likic
		"""
		
		return self.__time_list[:]
		
	@property
	def time_step(self):
		
		"""
		Returns the time step

		:return: Time step
		:rtype: float

		:author: Lewis Lee
		:author: Vladimir Likic
		"""
		
		return self.__time_step
	
	def __calc_time_step(self, time_list):
		
		"""
		Calculates the time step

		:param time_list: A list of retention times
		:type time_list: list

		:return: Time step value
		:rtype: float

		:author: Lewis Lee
		:author: Vladimir Likic
		"""
		
		td_list = []
		for ii in range(len(time_list) - 1):
			td = time_list[ii + 1] - time_list[ii]
			td_list.append(td)
		
		td_array = numpy.array(td_list)
		time_step = td_array.mean()
		
		return time_step
		
	def write(self, file_name, minutes=False, formatting=True):
		
		"""
		Writes the ion chromatogram to the specified file

		:param file_name: Output file name
		:type file_name: StringType
		:param minutes: A boolean value indicating whether to write
			time in minutes
		:type minutes: bool
		:param formatting: A boolean value indicating whether to
			format the numbers in the output (default True)
		:type minutes: bool

		:return: none
		:rtype: NoneType

		:author: Lewis Lee
		:author: Vladimir Likic
		"""
		
		if not isinstance(file_name, str):
			raise TypeError("'file_name' must be a string")
		
		fp = open_for_writing(file_name)
		
		time_list = copy.deepcopy(self.__time_list)
		
		if minutes:
			for ii in range(len(time_list)):
				time_list[ii] = time_list[ii] / 60.0
		
		for ii in range(len(time_list)):
			if formatting:
				fp.write("%8.4f %#.6e\n" % (time_list[ii], self.__ia[ii]))
			else:
				fp.write("{} {}\n".format(time_list[ii], self.__ia[ii]))
		
		close_for_writing(fp)

