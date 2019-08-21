"""
Class to model MS Scan data
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


import copy

import deprecation
from pyms import __version__

from pyms.Utils.Utils import is_list


class Scan(object):
	"""
	Generic object for a single Scan's raw data

	:author: Qiao Wang
	:author: Andrew Isaac
	:author: Vladimir Likic
	:author: Dominic Davis-Foster (type assertions and properties)
	"""
	
	def __init__(self, mass_list, intensity_list):
		
		"""
		Initialize the Scan data

		:param mass_list: mass values
		:type mass_list: list

		:param intensity_list: intensity values
		:type intensity_list: list

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""
		# print "mass_list[0]",mass_list[0]
		if not is_list(mass_list) or not isinstance(mass_list[0], (int, float)):
			raise TypeError("'mass_list' must be a list of numbers")
		if not is_list(intensity_list) or not isinstance(intensity_list[0], (int, float)):
			raise TypeError("'intensity_list' must be a list of numbers")
		
		self.__mass_list = mass_list
		self.__intensity_list = intensity_list
		self.__min_mass = min(mass_list)
		self.__max_mass = max(mass_list)
	
	def __len__(self):
		"""
		Returns the length of the Scan object

		:return: Length of Scan
		:rtype: int

		:author: Andrew Isaac
		"""
		
		return len(self.__mass_list)
	
	def __eq__(self, other):
		if isinstance(other, self.__class__):
			return self.intensity_list == other.intensity_list \
				   and self.mass_list == other.mass_list
		return NotImplemented
	
	def __copy__(self):
		return Scan(self.mass_list, self.intensity_list)
	
	def __deepcopy__(self, memodict={}):
		return self.__copy__()
	
	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'Scan.mass_list' instead")
	def get_mass_list(self):
		"""
		Returns the masses for the current scan

		:return: the masses
		:rtype: list

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""
		
		return self.mass_list
	
	@property
	def mass_list(self):
		"""
		Returns the masses for the current scan

		:return: the masses
		:rtype: list

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""
		
		return self.__mass_list[:]
	
	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'Scan.intensity_list' instead")
	def get_intensity_list(self):
		"""
		Returns the intensities for the current scan

		:return: the intensities
		:rtype: list

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""
		
		return self.intensity_list
	
	@property
	def intensity_list(self):
		"""
		Returns the intensities for the current scan

		:return: the intensities
		:rtype: list

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""
		
		return self.__intensity_list[:]
	
	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'Scan.min_mass' instead")
	def get_min_mass(self):
		
		"""
		Returns the minimum m/z value in the scan

		:return: Minimum m/z
		:rtype: Float

		:author: Andrew Isaac
		"""
		
		return self.min_mass
	
	@property
	def min_mass(self):
		
		"""
		Returns the minimum m/z value in the scan

		:return: Minimum m/z
		:rtype: Float

		:author: Andrew Isaac
		"""
		
		return self.__min_mass
	
	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'Scan.max_mass' instead")
	def get_max_mass(self):
		
		"""
		Returns the maximum m/z value in the scan

		:return: Maximum m/z
		:rtype: Float

		:author: Andrew Isaac
		"""
		
		return self.max_mass
	
	@property
	def max_mass(self):
		
		"""
		Returns the maximum m/z value in the scan

		:return: Maximum m/z
		:rtype: Float

		:author: Andrew Isaac
		"""
		
		return self.__max_mass

