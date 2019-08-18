"""
Class to model MS Scan data
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

import deprecation
from pyms import __version__

from pyms.Utils.Utils import is_list


class Scan(object):
	"""
	:summary: Generic object for a single Scan's raw data

	:author: Qiao Wang
	:author: Andrew Isaac
	:author: Vladimir Likic
	:author: Dominic Davis-Foster (type assertions and properties)
	"""
	
	def __init__(self, mass_list, intensity_list):
		
		"""
		:summary: Initialize the Scan data

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
		:summary: Returns the length of the Scan object

		:return: Length of Scan
		:rtype: int

		:author: Andrew Isaac
		"""
		
		return len(self.__mass_list)
	
	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'Scan.mass_list' instead")
	def get_mass_list(self):
		"""
		:summary: Returns the masses for the current scan

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
		:summary: Returns the masses for the current scan

		:return: the masses
		:rtype: list

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""
		
		return copy.deepcopy(self.__mass_list)
	
	# def __get_mass_list(self):
	#	return self.__mass_list
	
	# mass_list = property(__get_mass_list)
	
	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'Scan.intensity_list' instead")
	def get_intensity_list(self):
		"""
		:summary: Returns the intensities for the current scan

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
		:summary: Returns the intensities for the current scan

		:return: the intensities
		:rtype: list

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""
		
		return copy.deepcopy(self.__intensity_list)
	
	# def __get_intensity_list(self):
	#	return self.__intensity_list
	
	# intensity_list = property(__get_intensity_list)
	
	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'Scan.min_mass' instead")
	def get_min_mass(self):
		
		"""
		:summary: Returns the minimum m/z value in the scan

		:return: Minimum m/z
		:rtype: Float

		:author: Andrew Isaac
		"""
		
		return self.min_mass
	
	@property
	def min_mass(self):
		
		"""
		:summary: Returns the minimum m/z value in the scan

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
		:summary: Returns the maximum m/z value in the scan

		:return: Maximum m/z
		:rtype: Float

		:author: Andrew Isaac
		"""
		
		return self.max_mass
	
	@property
	def max_mass(self):
		
		"""
		:summary: Returns the maximum m/z value in the scan

		:return: Maximum m/z
		:rtype: Float

		:author: Andrew Isaac
		"""
		
		return self.__max_mass

