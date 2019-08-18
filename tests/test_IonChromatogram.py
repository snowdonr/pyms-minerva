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
import copy

import numpy

from pyms.GCMS.Class import IonChromatogram


def test_IonChromatogram(im, tic):
	# get the first ion chromatogram of the IntensityMatrix
	ic = im.get_ic_at_index(0)
	assert isinstance(ic, IonChromatogram)
	assert not ic.is_tic()
	
	# get the ion chromatogram for m/z = 73
	ic = im.get_ic_at_mass(73)
	assert isinstance(ic, IonChromatogram)
	assert not ic.is_tic()
	
	assert isinstance(tic, IonChromatogram)
	assert tic.is_tic()
	
	# Errors
	for type in [test_string, test_int, test_float, test_list_strs, test_list_ints, test_dict, test_tuple]:
		with pytest.raises(TypeError):
			IonChromatogram(type, tic.time_list)
	for type in [test_string, test_int, test_float, test_list_strs, test_dict]:
		with pytest.raises(TypeError):
			IonChromatogram(tic.intensity_array, type)
	for type in [test_string, test_list_strs, test_list_ints, test_dict, test_tuple]:
		with pytest.raises(TypeError):
			IonChromatogram(tic.intensity_array, tic.time_list, mass=type)
	
	with pytest.raises(ValueError):
		IonChromatogram(tic.intensity_array, test_list_ints)


def test_len(tic):
	assert len(tic) == 2103


def test_subtract_ic(im):
	ic1 = im.get_ic_at_index(0)
	assert isinstance(ic1, IonChromatogram)
	
	ic2 = im.get_ic_at_index(1)
	
	ic3 = ic1 - ic2
	assert isinstance(ic3, IonChromatogram)


def test_equality(tic, im):
	assert tic == IonChromatogram(tic.intensity_array, tic.time_list)
	assert tic != im.get_ic_at_index(0)
	assert tic != test_string
	assert tic != test_int
	assert tic != test_float
	assert tic != test_list_ints
	assert tic != test_list_strs
	assert tic != test_dict
	assert tic != test_tuple
	

def test_get_intensity_at_index(tic):
	assert isinstance(tic.get_intensity_at_index(test_int), float)
	assert tic.get_intensity_at_index(test_int) == 421170.0
	
	# Errors
	for type in [test_string, test_float, test_list_strs, test_list_ints, test_dict, test_tuple]:
		with pytest.raises(TypeError):
			tic.get_intensity_at_index(type)
	
	with pytest.raises(IndexError):
		tic.get_intensity_at_index(-1)
	with pytest.raises(IndexError):
		tic.get_intensity_at_index(10000000)


def test_get_intensity_array(tic):
	with pytest.warns(DeprecationWarning):
		tic.get_intensity_array()


def test_intensity_array(tic, im):
	tic = copy.deepcopy(tic)

	assert isinstance(tic.intensity_array, numpy.ndarray)
	assert all(numpy.equal(IonChromatogram(tic.intensity_array, tic.time_list).intensity_array, tic.intensity_array))
	
	ic = im.get_ic_at_index(0)
	tic.intensity_array = ic.intensity_array
	assert all(numpy.equal(tic.intensity_array, ic.intensity_array))


def test_set_intensity_array(tic):
	tic = copy.deepcopy(tic)
	with pytest.warns(DeprecationWarning):
		tic.set_intensity_array(tic.intensity_array)


def test_get_time_at_index(tic):
	assert isinstance(tic.get_time_at_index(test_int), float)
	assert tic.get_time_at_index(test_int) == 1304.15599823
	
	# Errors
	with pytest.raises(TypeError):
		tic.get_time_at_index(test_string)
	with pytest.raises(TypeError):
		tic.get_time_at_index(12.34)
	with pytest.raises(TypeError):
		tic.get_time_at_index([1, 2, 3, 4])
	with pytest.raises(TypeError):
		tic.get_time_at_index({"a": 1, "b": 2, "c": 3, "d": 4})
	# tic.get_time_at_index(0)
	
	with pytest.raises(IndexError):
		tic.get_time_at_index(-1)
	with pytest.raises(IndexError):
		tic.get_time_at_index(10000000)


def test_get_time_list(tic):
	with pytest.warns(DeprecationWarning):
		tic.get_time_list()


def test_time_list(tic):
	assert isinstance(tic.time_list, list)
	assert isinstance(tic.time_list[0], float)
	assert tic.time_list[0] == 1.05200003833


def test_get_mass(im):
	ic = im.get_ic_at_index(0)
	with pytest.warns(DeprecationWarning):
		ic.get_mass()


def test_mass(tic, im):
	with pytest.warns(Warning):
		tic.mass
	
	ic = im.get_ic_at_index(0)
	assert isinstance(ic.mass, (int, float))
	assert ic.mass == 50.2516


def test_get_time_step(tic):
	with pytest.warns(DeprecationWarning):
		tic.get_time_step()


def test_time_step(tic):
	assert isinstance(tic.time_step, float)
	assert tic.time_step == 1.0560000035830972
	

def test_get_index_at_time(tic):
	assert isinstance(tic.get_index_at_time(12), int)
	assert tic.get_index_at_time(12) == 10
	
	# Errors
	for type in [test_string, test_list_ints, test_list_strs, test_dict, test_tuple]:
		with pytest.raises(TypeError):
			tic.get_index_at_time(type)
	
	with pytest.raises(IndexError):
		tic.get_index_at_time(-1)
	
	with pytest.raises(IndexError):
		tic.get_index_at_time(1000000)


def test_write(tic):
	tic.write("output/tic.dat",minutes=False, formatting=False)
	
	from pyms.Utils.IO import open_for_reading, close_for_reading
	
	fp = open_for_reading("output/tic.dat")
	
	for line, ii in zip(fp.readlines(), range(len(tic.time_list))):
		assert line == "{} {}\n".format(tic.time_list[ii], tic.intensity_array[ii])
		
	close_for_reading(fp)
	
	tic.write("output/tic_minutes.dat",minutes=True, formatting=False)
	
	fp = open_for_reading("output/tic_minutes.dat")
	
	for line, ii in zip(fp.readlines(), range(len(tic.time_list))):
		assert line == "{} {}\n".format(tic.time_list[ii]/60.0, tic.intensity_array[ii])
		
	close_for_reading(fp)
	
	tic.write("output/tic_formatting.dat", minutes=False)
	
	fp = open_for_reading("output/tic_formatting.dat")
	
	for line, ii in zip(fp.readlines(), range(len(tic.time_list))):
		assert line == "%8.4f %#.6e\n" % (tic.time_list[ii], tic.intensity_array[ii])
		
	close_for_reading(fp)
	
	for type in [test_dict, test_list_strs, test_list_ints, test_tuple, test_int, test_float]:
		with pytest.raises(TypeError):
			tic.write(type)

