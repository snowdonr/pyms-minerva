#############################################################################
#                                                                           #
#    PyMS software for processing of metabolomic mass-spectrometry data     #
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

import pytest
from tests.constants import *

from copy import deepcopy

from pyms.GCMS.Class import GCMS_data
from pyms.IonChromatogram import IonChromatogram
from pyms.Scan import Scan


def test_GCMS_data(data):
	assert isinstance(data, GCMS_data)
	assert len(data) == 2103
	
	GCMS_data(data.time_list, data.scan_list)
	
	# Errors
	for type in [test_string, test_float, test_int, test_list_strs, test_dict]:
		with pytest.raises(TypeError):
			GCMS_data(type, data.scan_list)
			
	for type in [test_string, test_float, test_int, test_list_ints, test_list_strs, test_dict, test_tuple]:
		with pytest.raises(TypeError):
			GCMS_data(data.time_list, type)
		

def test_get_min_mass(data):
	with pytest.warns(DeprecationWarning):
		data.get_min_mass()


def test_min_mass(data):
	assert isinstance(data.min_mass, float)
	# minimum mass found in all data
	assert data.min_mass == 50.2516


def test_get_max_mass(data):
	with pytest.warns(DeprecationWarning):
		data.get_max_mass()


def test_max_mass(data):
	# maximum mass found in all data
	assert isinstance(data.max_mass, float)
	assert data.max_mass == 499.6226


def test_get_index_at_time(data):
	# index of 400sec in time_list
	assert isinstance(data.get_index_at_time(400.0), int)
	assert data.get_index_at_time(400.0) == 378


def test_get_time_list(data):
	with pytest.warns(DeprecationWarning):
		data.get_time_list()


def test_time_list(data):
	time = data.time_list
	assert isinstance(time, list)
	# number of retention times
	assert len(time) == 2103
	# retention time of 1st scan:
	assert isinstance(time[0], float)
	assert time[0] == 1.05200003833


def test_get_scan_list(data):
	with pytest.warns(DeprecationWarning):
		data.get_scan_list()


def test_scan_list(data):
	# raw scans
	scans = data.scan_list
	
	assert isinstance(scans, list)
	assert isinstance(scans[0], Scan)
	assert len(scans[0]) == 101
	assert isinstance(scans[0].mass_list, list)
	# 1st mass value for 1st scan
	assert isinstance(scans[0].mass_list[0], float)
	assert scans[0].mass_list[0] == 52.0131
	
	assert isinstance(scans[0].intensity_list, list)
	# 1st intensity value for 1st scan
	assert isinstance(scans[0].intensity_list[0], float)
	assert scans[0].intensity_list[0] == 5356.0
	
	# minimum mass found in 1st scan
	assert isinstance(scans[0].min_mass, float)
	assert scans[0].min_mass == 52.0131
	
	# maximum mass found in 1st scan
	assert isinstance(scans[0].max_mass, float)
	assert scans[0].max_mass == 477.6667


def test_get_tic(data):
	with pytest.warns(DeprecationWarning):
		data.get_tic()


def test_tic(data):
	tic = data.tic
	assert isinstance(tic, IonChromatogram)
	# number of scans in TIC
	assert len(tic) == 2103
	assert len(tic) == len(data.time_list)
	
	# start time of TIC
	assert isinstance(tic.get_time_at_index(0), float)
	assert tic.get_time_at_index(0) == 1.05200003833
	assert isinstance(tic.get_index_at_time(1.05200003833), int)
	assert tic.get_index_at_time(1.05200003833) == 0
	assert tic.get_index_at_time(2) == 1
	
	assert isinstance(tic.get_intensity_at_index(44), float)
	assert tic.get_intensity_at_index(44) == 3113615.0

	assert isinstance(tic.time_list, list)
	assert tic.time_list[0] == 1.05200003833
	
	assert isinstance(tic.time_step, float)
	
	assert isinstance(tic.is_tic(), bool)
	assert tic.is_tic()


def test_trim(data):
	# time
	trimmed = deepcopy(data)
	trimmed.trim("6.5m", "21m")
	
	assert trimmed.min_mass == 50.2516
	assert trimmed.max_mass == 499.6226
	
	time = trimmed.time_list
	assert len(time) == 825
	assert time[0] == 390.715999603
	assert trimmed.get_index_at_time(400.0) == 9
	
	scans = trimmed.scan_list
	assert scans[0].mass_list[0] == 51.0066

	# Scans
	trimmed = deepcopy(data)
	trimmed.trim(1000, 2000)
	
	assert trimmed.min_mass == 50.2516
	assert trimmed.max_mass == 499.6226
	
	time = trimmed.time_list
	assert len(time) == 1002
	assert time[0] == 1055.99601746
	assert trimmed.get_index_at_time(1500.0) == 420
	
	scans = trimmed.scan_list
	assert scans[0].mass_list[0] == 50.8808


def test_info(capsys, data):
	data.info()
	captured = capsys.readouterr()
	assert captured.out == """ Data retention time range: 0.018 min -- 37.013 min
 Time step: 1.056 s (std=0.000 s)
 Number of scans: 2103
 Minimum m/z measured: 50.252
 Maximum m/z measured: 499.623
 Mean number of m/z values per scan: 98
 Median number of m/z values per scan: 98
"""


def test_write(data):
	raise NotImplementedError


def test_write_intensities_stream(data):
	raise NotImplementedError
