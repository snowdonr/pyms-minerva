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

from pyms.GCMS.Class import MassSpectrum


def test_MassSpectrum(ms):
	assert isinstance(ms, MassSpectrum)
	
	assert isinstance(ms.mass_list, list)
	assert isinstance(ms.mass_spec, list)
	
	# Errors
	for type in [test_string, test_int, test_list_strs, test_dict]:
		with pytest.raises(TypeError):
			MassSpectrum(type, ms.mass_spec)
	
	for type in [test_string, test_int, test_list_strs, test_dict]:
		with pytest.raises(TypeError):
			MassSpectrum(ms.mass_list, type)
	
	with pytest.raises(ValueError):
		MassSpectrum(ms.mass_list, test_list_ints)


def test_len(ms):
	assert len(ms) == 450


def test_equality(im, ms):
	assert ms != im.get_ms_at_index(1234)
	assert ms == MassSpectrum(ms.mass_list, ms.mass_spec)
	assert ms != test_list_ints
	assert ms != test_list_strs
	assert ms != test_tuple
	assert ms != test_string
	assert ms != test_int
	assert ms != test_float

