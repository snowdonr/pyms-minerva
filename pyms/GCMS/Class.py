"""
Classes to model GC-MS data
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


import numpy
import math
import copy

import deprecation
from pyms import __version__

from pyms.Utils.Error import pymsError
from pyms.Utils.Utils import is_str, is_list
from pyms.Utils.IO import open_for_writing, close_for_writing, save_data
from pyms.Utils.Math import mean, std, median
from pyms.Utils.Time import time_str_secs
from pyms.Scan import Scan
from pyms.IonChromatogram import IonChromatogram
from pyms.IntensityMatrix import IntensityMatrix
from pyms.MassSpectrum import MassSpectrum


class GCMS_data(object):
	"""
	:summary: Generic object for GC-MS data. Contains raw data
		as a list of scans and times

	:author: Qiao Wang
	:author: Andrew Isaac
	:author: Vladimir Likic
	:author: Dominic Davis-Foster (type assertions and properties)
	"""

	def __init__(self, time_list, scan_list):
		"""
		:summary: Initialize the GC-MS data

		:param time_list: List of scan retention times
		:type time_list: list
		:param scan_list: List of Scan objects
		:type scan_list: list

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""

		if not is_list(time_list) or not isinstance(time_list[0], (int, float)):
			raise TypeError("'time_list' must be a list of numbers")

		if not is_list(scan_list) or not isinstance(scan_list[0], Scan):
			raise TypeError("'scan_list' must be a list of Scan objects")

		self.__set_time(time_list)
		self.__scan_list = scan_list
		self.__set_min_max_mass()
		self.__calc_tic()

	def __len__(self):
		"""
		:summary: Returns the length of the data object,
			defined as the number of scans

		:return: Number of scans
		:rtype: int

		:author: Vladimir Likic
		"""

		return len(self.__scan_list)

	def __set_time(self, time_list):
		"""
		:summary: Sets time-related properties of the data

		:param time_list: List of retention times
		:type time_list: list

		:author: Vladimir Likic
		"""

		# calculate the time step, its spreak, and along the way
		# check that retention times are increasing
		time_diff_list = []

		for ii in range(len(time_list)-1):
			t1 = time_list[ii]
			t2 = time_list[ii+1]
			if not t2 > t1:
				raise pymsError("problem with retention times detected")
			time_diff = t2 - t1
			time_diff_list.append(time_diff)

		time_step = mean(time_diff_list)
		time_step_std = std(time_diff_list)

		self.__time_list = time_list
		self.__time_step = time_step
		self.__time_step_std = time_step_std
		self.__min_rt = min(time_list)
		self.__max_rt = max(time_list)

	def __set_min_max_mass(self):
		"""
		:summary: Sets the min and max mass value

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""

		mini = self.__scan_list[0].min_mass
		maxi = self.__scan_list[0].max_mass
		for scan in self.__scan_list:
			tmp_mini = scan.min_mass
			tmp_maxi = scan.max_mass
			if tmp_mini < mini:
				mini = tmp_mini
			if tmp_maxi > maxi:
				maxi = tmp_maxi
		self.__min_mass = mini
		self.__max_mass = maxi

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'GCMS>min_mass' instead")
	def get_min_mass(self):
		"""
		:summary: Get the min mass value over all scans

		:return: The minimum mass of all the data
		:rtype: float

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""

		return self.min_mass
	
	@property
	def min_mass(self):
		"""
		:summary: Get the min mass value over all scans

		:return: The minimum mass of all the data
		:rtype: float

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""

		return self.__min_mass

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'GCMS>max_mass' instead")
	def get_max_mass(self):
		"""
		:summary: Get the max mass value over all scans

		:return: The maximum mass of all the data
		:rtype: float

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""

		return self.max_mass
	
	@property
	def max_mass(self):
		"""
		:summary: Get the max mass value over all scans

		:return: The maximum mass of all the data
		:rtype: float

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""

		return self.__max_mass

	def get_index_at_time(self, time):

		"""
		:summary: Returns the nearest index corresponding to the given
			time

		:param time: Time in seconds
		:type time: float

		:return: Nearest index corresponding to given time
		:rtype: int

		:author: Lewis Lee
		:author: Tim Erwin
		:author: Vladimir Likic
		"""

		if not isinstance(time, (int, float)):
			raise TypeError("'time' must be a number")

		if (time < self.__min_rt) or (time > self.__max_rt):
			raise IndexError(f"time {time:.2f} is out of bounds (min: {self.__min_rt:.2f}, max: {self.__max_rt:.2f})")

		time_list = self.__time_list
		time_diff_min = self.__max_rt
		ix_match = None

		for ix in range(len(time_list)):

			time_diff = math.fabs(time-time_list[ix])

			if time_diff < time_diff_min:
				ix_match = ix
				time_diff_min = time_diff

		return ix_match
	
	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'GCMS.time_list' instead")
	def get_time_list(self):
		"""
		:summary: Returns the list of each scan retention time

		:return: A list of each scan retention time
		:rtype: list

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""

		return self.time_list

	@property
	def time_list(self):
		"""
		:summary: Returns the list of each scan retention time

		:return: A list of each scan retention time
		:rtype: list

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""

		return copy.deepcopy(self.__time_list)

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'GCMS.scan_list' instead")
	def get_scan_list(self):
		"""
		:summary: Return a list of the scan objects

		:return: A list of scan objects
		:rtype: list

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""

		return self.scan_list
	
	@property
	def scan_list(self):
		"""
		:summary: Return a list of the scan objects

		:return: A list of scan objects
		:rtype: list

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""

		return copy.deepcopy(self.__scan_list)

	#def __get_scan_list(self):
	#	return self.__scan_list

	#scan_list = property(__get_scan_list)


	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'GCMS.tic' instead")
	def get_tic(self):

		"""
		:summary: Returns the total ion chromatogram

		:return: Total ion chromatogram
		:rtype: pyms.GCMS.Class.IonChromatogram

		:author: Andrew Isaac
		"""

		return self.tic
	
	@property
	def tic(self):

		"""
		:summary: Returns the total ion chromatogram

		:return: Total ion chromatogram
		:rtype: pyms.GCMS.Class.IonChromatogram

		:author: Andrew Isaac
		"""

		return self.__tic
	
	def __calc_tic(self):
		"""
		:summary: Calculate the total ion chromatogram

		:return: Total ion chromatogram
		:rtype: pyms.GCMS.Class.IonChromatogram

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""

		intensities = []
		for scan in self.__scan_list:
			intensities.append(sum(scan.intensity_list))
		ia = numpy.array(intensities)
		rt = copy.deepcopy(self.__time_list)
		tic = IonChromatogram(ia, rt)

		self.__tic = tic

	def trim(self, begin=None, end=None):
		"""
		:summary: trims data in the time domain

		:param begin: begin parameter designating start time or
			scan number
		:type begin: int or StrType
		:param end: end parameter designating start time or
			scan number
		:type end: int or StrType

			The arguments 'begin' and 'end' can be either integers
			(in which case they are taken as the first/last scan
			number for trimming) or strings in which case they are
			treated as time strings and converted to scan numbers.

			At present both 'begin' and 'end' must be of the same
			type, either both scan numbers or time strings.

		:author: Vladimir Likic
		"""

		# trim called with defaults, or silly arguments
		if begin == None and end == None:
			print("Nothing to do.")
			return # exit immediately

		N = len(self.__scan_list)

		# process 'begin' and 'end'
		if begin == None:
			first_scan = 0
		elif isinstance(begin, (int, float)):
			first_scan = begin-1
		elif is_str(begin):
			time = time_str_secs(begin)
			first_scan = self.get_index_at_time(time) + 1
		else:
			raise ValueError("invalid 'begin' argument")

		if end == None:
			last_scan = N-1
		elif isinstance(end, (int, float)):
			last_scan = end
		elif is_str(end):
			time = time_str_secs(end)
			last_scan = self.get_index_at_time(time) + 1
		else:
			raise ValueError("invalid 'end' argument")

		# sanity checks
		if not last_scan > first_scan:
			raise ValueError("last scan=%d, first scan=%d" % (last_scan, first_scan))
		elif first_scan < 0:
			raise ValueError("scan number must be greater than one")
		elif last_scan > N-1:
			raise ValueError("last scan=%d, total number of scans=%d" % (last_scan, N))

		print("Trimming data to between %d and %d scans" % \
				(first_scan+1, last_scan+1))

		scan_list_new = []
		time_list_new = []
		for ii in range(len(self.__scan_list)):
			if ii >= first_scan and ii <= last_scan:
				scan = self.__scan_list[ii]
				time = self.__time_list[ii]
				scan_list_new.append(scan)
				time_list_new.append(time)

		# update info
		self.__scan_list = scan_list_new
		self.__set_time(time_list_new)
		self.__set_min_max_mass()
		self.__calc_tic()

	def info(self, print_scan_n=False):
		"""
		:summary: Prints some information about the data

		:param print_scan_n: If set to True will print the number
			of m/z values in each scan
		:type print_scan_n: BooleanType

		:author: Vladimir Likic
		"""

		# print the summary of simply attributes
		print(" Data retention time range: %.3f min -- %.3f min" % (self.__min_rt/60.0, self.__max_rt/60))
		print(" Time step: %.3f s (std=%.3f s)" % ( self.__time_step, self.__time_step_std ))
		print(" Number of scans: %d" % ( len(self.__scan_list) ))
		print(" Minimum m/z measured: %.3f" % ( self.__min_mass ))
		print(" Maximum m/z measured: %.3f" % ( self.__max_mass ))

		# calculate median number of m/z values measured per scan
		n_list = []
		for ii in range(len(self.__scan_list)):
			scan = self.__scan_list[ii]
			n = len(scan)
			n_list.append(n)
			if print_scan_n: print(n)
		mz_mean = mean(n_list)
		mz_median = median(n_list)
		print(" Mean number of m/z values per scan: %d" % ( mz_mean ))
		print(" Median number of m/z values per scan: %d" % ( mz_median ))

	def write(self, file_root):

		"""
		:summary: Writes the entire raw data to two files, one
			'file_root'.I.csv (intensities) and 'file_root'.mz.csv
			(m/z values).

			This method writes two CSV files, containing intensities
			and corresponding m/z values. In general these are not
			two-dimensional matrices, because different scans may
			have different number of m/z values recorded.

		:param file_root: The root for the output file names
		:type file_root: StringType

		:author: Vladimir Likic
		"""

		if not is_str(file_root):
			raise TypeError("'file_root' must be a string")

		file_name1 = file_root + ".I.csv"
		file_name2 = file_root + ".mz.csv"

		print(" -> Writing intensities to '%s'" % ( file_name1 ))
		print(" -> Writing m/z values to '%s'" % ( file_name2 ))

		fp1 = open_for_writing(file_name1)
		fp2 = open_for_writing(file_name2)

		for ii in range(len(self.__scan_list)):

			scan = self.__scan_list[ii]

			intensity_list = scan.get_intensity_list()
			mass_list = scan.get_mass_list()

			for ii in range(len(intensity_list)):
				v = intensity_list[ii]
				if ii == 0:
					fp1.write("%.4f" % (v))
				else:
					fp1.write(",%.4f" % (v))
			fp1.write("\n")

			for ii in range(len(mass_list)):
				v = mass_list[ii]
				if ii == 0:
					fp2.write("%.4f" % (v))
				else:
					fp2.write(",%.4f" % (v))
			fp2.write("\n")

		close_for_writing(fp1)
		close_for_writing(fp2)

	def write_intensities_stream(self, file_name):

		"""
		:summary: Writes all intensities to a file

		:param file_name: Output file name
		:type file_name: StringType

		This function loop over all scans, and for each scan
		writes intensities to the file, one intenisity per
		line. Intensities from different scans are joined
		without any delimiters.

		:author: Vladimir Likic
		"""

		if not is_str(file_name):
			raise TypeError("'file_name' must be a string")

		N = len(self.__scan_list)

		print(" -> Writing scans to a file")
		
		fp = open_for_writing(file_name)

		for ii in range(len(self.__scan_list)):
			scan = self.__scan_list[ii]
			intensities = scan.get_intensity_list()
			for I in intensities:
				fp.write(f"{I:8.4f}\n")

		close_for_writing(fp)

## get_ms_at_time()
