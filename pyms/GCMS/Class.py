"""
Class to model GC-MS data
"""

################################################################################
#                                                                              #
#    PyMassSpec software for processing of metabolomic mass-spectrometry data  #
#    Copyright (C) 2005-2012 Vladimir Likic                                    #
#    Copyright (C) 2019 Dominic Davis-Foster                                   #
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


import copy
import math
import pathlib

import numpy
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
	Generic object for GC-MS data. Contains raw data
		as a list of scans and times

	:param time_list: List of scan retention times
	:type time_list: list
	:param scan_list: List of Scan objects
	:type scan_list: list

	:author: Qiao Wang
	:author: Andrew Isaac
	:author: Vladimir Likic
	:author: Dominic Davis-Foster (type assertions and properties)
	"""
	
	def __init__(self, time_list, scan_list):
		"""
		Initialize the GC-MS data
		"""
		
		if not is_list(time_list) or not isinstance(time_list[0], (int, float)):
			raise TypeError("'time_list' must be a list of numbers")
		
		if not is_list(scan_list) or not isinstance(scan_list[0], Scan):
			raise TypeError("'scan_list' must be a list of Scan objects")
		
		self.__time_list = time_list
		self.__scan_list = scan_list
		self.__set_time()
		self.__set_min_max_mass()
		self.__calc_tic()
	
	def __eq__(self, other):
		# 		self.__time_step = time_step
		# 		self.__time_step_std = time_step_std
		
		if isinstance(other, self.__class__):
			return self.scan_list == other.scan_list \
				   and self.time_list == other.time_list
		return NotImplemented
	
	def __len__(self):
		"""
		Returns the length of the data object, defined as the number of scans

		:return: Number of scans
		:rtype: int

		:author: Vladimir Likic
		"""
		
		return len(self.__scan_list)
	
	def __calc_tic(self):
		"""
		Calculate the total ion chromatogram

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
	
	def __set_time(self):
		"""
		Sets time-related properties of the data

		:author: Vladimir Likic
		"""
		
		# calculate the time step, its spreak, and along the way
		# check that retention times are increasing
		time_diff_list = []
		
		for index, t1 in enumerate(self.__time_list):
			t2 = self.__time_list[index + 1]
			if not t2 > t1:
				raise pymsError("problem with retention times detected")
			time_diff = t2 - t1
			time_diff_list.append(time_diff)
		
		time_step = mean(time_diff_list)
		time_step_std = std(time_diff_list)
		
		self.__time_step = time_step
		self.__time_step_std = time_step_std
		self.__min_rt = min(self.__time_list)
		self.__max_rt = max(self.__time_list)
	
	def __set_min_max_mass(self):
		"""
		Sets the min and max mass value

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
		
	
	
	def get_index_at_time(self, time):
		"""
		Returns the nearest index corresponding to the given time

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
			
			time_diff = math.fabs(time - time_list[ix])
			
			if time_diff < time_diff_min:
				ix_match = ix
				time_diff_min = time_diff
		
		return ix_match
	
	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'GCMS_data.max_mass' instead")
	def get_max_mass(self):
		"""
		Get the max mass value over all scans

		.. deprecated:: 2.1.2
			Use :attr:`pyms.GCMS.Class.GCMS_data.max_mass` instead.

		:return: The maximum mass of all the data
		:rtype: float

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""
		
		return self.max_mass
	
	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'GCMS_data.min_mass' instead")
	def get_min_mass(self):
		"""
		Get the min mass value over all scans
		
		.. deprecated:: 2.1.2
			Use :attr:`pyms.GCMS.Class.GCMS_data.min_mass` instead.

		:return: The minimum mass of all the data
		:rtype: float

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""
		
		return self.min_mass
	
	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'GCMS.scan_list' instead")
	def get_scan_list(self):
		"""
		Return a list of the scan objects

		.. deprecated:: 2.1.2
			Use :attr:`pyms.GCMS.Class.GCMS_data.scan_list` instead.

		:return: A list of scan objects
		:rtype: list

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""
		
		return self.scan_list
	
	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'GCMS.tic' instead")
	def get_tic(self):
		"""
		Returns the total ion chromatogram

		.. deprecated:: 2.1.2
			Use :attr:`pyms.GCMS.Class.GCMS_data.tic` instead.

		:return: Total ion chromatogram
		:rtype: pyms.GCMS.Class.IonChromatogram

		:author: Andrew Isaac
		"""
		
		return self.tic

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'GCMS.time_list' instead")
	def get_time_list(self):
		"""
		Returns the list of each scan retention time

		.. deprecated:: 2.1.2
			Use :attr:`pyms.GCMS.Class.GCMS_data.time_list` instead.

		:return: A list of each scan retention time
		:rtype: list

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""
		
		return self.time_list
	
	def info(self, print_scan_n=False):
		"""
		Prints some information about the data

		:param print_scan_n: If set to True will print the number of m/z values in each scan
		:type print_scan_n: bool, optional

		:author: Vladimir Likic
		"""
		
		# print the summary of simply attributes
		print(f" Data retention time range: {self.__min_rt / 60.0:.3f} min -- {self.__max_rt / 60:.3f} min")
		print(f" Time step: {self.__time_step:.3f} s (std={self.__time_step_std:.3f} s)")
		print(f" Number of scans: {len(self.__scan_list):d}")
		print(f" Minimum m/z measured: {self.__min_mass:.3f}")
		print(f" Maximum m/z measured: {self.__max_mass:.3f}")
		
		# calculate median number of m/z values measured per scan
		n_list = []
		for ii in range(len(self.__scan_list)):
			scan = self.__scan_list[ii]
			n = len(scan)
			n_list.append(n)
			if print_scan_n: print(n)
		mz_mean = mean(n_list)
		mz_median = median(n_list)
		print(f" Mean number of m/z values per scan: {mz_mean:d}")
		print(f" Median number of m/z values per scan: {mz_median:d}")
	
	@property
	def max_mass(self):
		"""
		Get the max mass value over all scans

		:return: The maximum mass of all the data
		:rtype: float

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""
		
		return self.__max_mass

	@property
	def min_mass(self):
		"""
		Get the min mass value over all scans

		:return: The minimum mass of all the data
		:rtype: float

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""
		
		return self.__min_mass
	
	@property
	def scan_list(self):
		"""
		Return a list of the scan objects

		:return: A list of scan objects
		:rtype: list

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""
		
		return copy.deepcopy(self.__scan_list)
	
	@property
	def tic(self):
		"""
		Returns the total ion chromatogram

		:return: Total ion chromatogram
		:rtype: pyms.GCMS.Class.IonChromatogram

		:author: Andrew Isaac
		"""
		
		return self.__tic
	
	@property
	def time_list(self):
		"""
		Returns the list of each scan retention time

		:return: A list of each scan retention time
		:rtype: list

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""
		
		return self.__time_list[:]
	
	def trim(self, begin=None, end=None):
		"""
		trims data in the time domain
		
		The arguments ``begin`` and ``end`` can be either integers (in which case
		they are taken as the first/last scan number for trimming) or strings
		in which case they are treated as time strings and converted to scan
		numbers.

		At present both ``begin`` and ``end`` must be of the same type, either both
		scan numbers or time strings.
		
		At least one of ``begin`` and ``end`` is required

		:param begin: begin parameter designating start time or scan number
		:type begin: int or str, optional
		:param end: end parameter designating start time or scan number
		:type end: int or str, optional

		:author: Vladimir Likic
		"""
		
		# trim called with defaults, or silly arguments
		if begin is None and end is None:
			raise SyntaxError("At least one of 'begin' and 'end' is required")
		
		N = len(self.__scan_list)
		
		# process 'begin' and 'end'
		if begin is None:
			first_scan = 0
		elif isinstance(begin, (int, float)):
			first_scan = begin - 1
		elif is_str(begin):
			time = time_str_secs(begin)
			first_scan = self.get_index_at_time(time) + 1
		else:
			raise ValueError("invalid 'begin' argument")
		
		if end is None:
			last_scan = N - 1
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
		elif last_scan > N - 1:
			raise ValueError("last scan=%d, total number of scans=%d" % (last_scan, N))
		
		print("Trimming data to between %d and %d scans" % \
			  (first_scan + 1, last_scan + 1))
		
		scan_list_new = []
		time_list_new = []
		for ii in range(len(self.__scan_list)):
			if first_scan <= ii <= last_scan:
				scan = self.__scan_list[ii]
				time = self.__time_list[ii]
				scan_list_new.append(scan)
				time_list_new.append(time)
		
		# update info
		self.__scan_list = scan_list_new
		self.__time_list = time_list_new
		self.__set_time()
		self.__set_min_max_mass()
		self.__calc_tic()
		
	def write(self, file_root):
		"""
		Writes the entire raw data to two files, one
			'file_root'.I.csv (intensities) and 'file_root'.mz.csv
			(m/z values).

			This method writes two CSV files, containing intensities
			and corresponding m/z values. In general these are not
			two-dimensional matrices, because different scans may
			have different number of m/z values recorded.

		:param file_root: The root for the output file names
		:type file_root: str or pathlib.Path

		:author: Vladimir Likic
		:author: Dominic Davis-Foster (pathlib support)
		"""
		
		if not isinstance(file_root, (str, pathlib.Path)):
			raise TypeError("'file_root' must be a string or a pathlib.Path object")
		
		if not isinstance(file_root, pathlib.Path):
			file_root = pathlib.Path(file_root)
		
		if not file_root.parent.is_dir():
			file_root.parent.mkdir(parents=True)
		
		
		file_name1 = str(file_root) + ".I.csv"
		file_name2 = str(file_root) + ".mz.csv"
		
		print(f" -> Writing intensities to '{file_name1}'")
		print(f" -> Writing m/z values to '{file_name2}'")
		
		fp1 = open_for_writing(file_name1)
		fp2 = open_for_writing(file_name2)
		
		for scan in self.__scan_list:
			
			for index, intensity in enumerate(scan.intensity_list):
				if index == 0:
					fp1.write(f"{intensity:.4f}")
				else:
					fp1.write(f",{intensity:.4f}")
			fp1.write("\n")
			
			for index, mass in enumerate(scan.mass_list):
				if index == 0:
					fp2.write(f"{mass:.4f}")
				else:
					fp2.write(f",{mass:.4f}")
			fp2.write("\n")
		
		close_for_writing(fp1)
		close_for_writing(fp2)
	
	def write_intensities_stream(self, file_name):
		"""
		Writes all intensities to a file

		This function loop over all scans, and for each scan
		writes intensities to the file, one intensity per
		line. Intensities from different scans are joined
		without any delimiters.

		:param file_name: Output file name
		:type file_name: str or pathlib.Path

		:author: Vladimir Likic
		:author: Dominic Davis-Foster (pathlib support)
		"""
		
		if not isinstance(file_name, (str, pathlib.Path)):
			raise TypeError("'file_name' must be a string or a pathlib.Path object")
		
		if not isinstance(file_name, pathlib.Path):
			file_name = pathlib.Path(file_name)
		
		if not file_name.parent.is_dir():
			file_name.parent.mkdir(parents=True)
		
		N = len(self.__scan_list)
		
		print(" -> Writing scans to a file")
		
		fp = file_name.open(file_name)
		
		for scan in self.__scan_list:
			intensities = scan.get_intensity_list()
			for I in intensities:
				fp.write(f"{I:8.4f}\n")
		
		fp.close()

## get_ms_at_time()
