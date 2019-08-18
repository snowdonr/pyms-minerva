"""
Classes to model a Mass Spectrum
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


from pyms.Utils.Utils import is_list


class MassSpectrum(object):
	"""
	:summary: Models a binned mass spectrum

	:author: Andrew Isaac
	:author: Qiao Wang
	:author: Vladimir Likic
	"""
	
	def __init__(self, mass_list, intensity_list):
		"""
		:summary: Initialise the MassSpectrum

		:param mass_list: List of binned masses
		:type mass_list: list
		:param intensity_list: List of binned intensities
		:type intensity_list: list

		:author: Andrew Isaac
		:author: Qiao Wang
		:author: Vladimir Likic
		"""
		
		if not is_list(mass_list) or not isinstance(mass_list[0], (int, float)):
			raise TypeError("'mass_list' must be a list of numbers")
		if not is_list(intensity_list) or \
				not isinstance(intensity_list[0], (int, float)):
			raise TypeError("'intensity_list' must be a list of numbers")
		if not len(mass_list) == len(intensity_list):
			raise ValueError("'mass_list' is not the same size as 'intensity_list'")
		
		# TODO: should these be public, or accessed through methods???
		self.mass_list = mass_list
		self.mass_spec = intensity_list
	
	def __len__(self):
		
		"""
		:summary: Length of the MassSpectrum

		:return: Length of the MassSpectrum (Number of bins)
		:rtype: int

		:author: Andrew Isaac
		:author: Qiao Wang
		:author: Vladimir Likic
		"""
		
		return len(self.mass_list)
	
	def __eq__(self, other):
		if isinstance(other, self.__class__):
			return self.mass_list == other.mass_list and self.mass_spec == other.mass_spec
		return NotImplemented

	def __copy__(self):
		return MassSpectrum(mass_list=self.mass_list[:], intensity_list=self.mass_spec[:])
	
	def __deepcopy__(self, memodict={}):
		return self.__copy__()
