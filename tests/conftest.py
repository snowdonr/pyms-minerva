import pytest

import os
import shutil
import sys

from pathlib import Path

from pyms.TopHat import tophat
from pyms.BillerBiemann import BillerBiemann, rel_threshold, num_ions_threshold
from pyms.GCMS.Function import build_intensity_matrix, build_intensity_matrix_i
from pyms.GCMS.IO.JCAMP.Function import JCAMP_reader
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.Peak.Function import peak_sum_area, peak_top_ion_areas
from pyms.Peak.Class import Peak
from pyms.Experiment.Class import Experiment

from copy import deepcopy
# TODO: consider using pytest-datadir


# From pytest-datadir (https://github.com/gabrielcnr/pytest-datadir)
def _win32_longpath(path):
	'''
	Helper function to add the long path prefix for Windows, so that shutil.copytree won't fail
	while working with paths with 255+ chars.
	'''
	if sys.platform == 'win32':
		# The use of os.path.normpath here is necessary since "the "\\?\" prefix to a path string
		# tells the Windows APIs to disable all string parsing and to send the string that follows
		# it straight to the file system".
		# (See https://docs.microsoft.com/pt-br/windows/desktop/FileIO/naming-a-file)
		return '\\\\?\\' + os.path.normpath(path)
	else:
		return path


@pytest.fixture("session")
def datadir():
	return Path(os.path.split(__file__)[0]) / "data"


@pytest.fixture(scope="session")
def data(datadir):
	print("data")
	#return JCAMP_reader(datadir / "ELEY_1_SUBTRACT.JDX")
	return JCAMP_reader(os.path.join("tests", "data", "ELEY_1_SUBTRACT.JDX"))


@pytest.fixture(scope="session")
def im(data):
	# build an intensity matrix object from the data
	return build_intensity_matrix(data)


@pytest.fixture(scope="session")
def tic(data):
	# get the TIC
	return deepcopy(data.tic)


@pytest.fixture(scope="session")
def im_i(data):
	# build an intensity matrix object from the data
	return build_intensity_matrix_i(data)


@pytest.fixture(scope="function")
def peak_list(im_i):
	im_i = deepcopy(im_i)
	
	# Intensity matrix size (scans, masses)
	n_scan, n_mz = im_i.size
	
	# noise filter and baseline correct
	for ii in range(n_mz):
		ic = im_i.get_ic_at_index(ii)
		ic_smooth = savitzky_golay(ic)
		ic_bc = tophat(ic_smooth, struct="1.5m")
		im_i.set_ic_at_index(ii, ic_bc)
	
	# Use Biller and Biemann technique to find apexing ions at a scan
	# default is maxima over three scans and not to combine with any neighbouring
	# scan.
	peak_list = BillerBiemann(im_i, points=9, scans=2)
	return peak_list
	

@pytest.fixture(scope="function")
def filtered_peak_list(im_i, peak_list):
	#peak_list = deepcopy(peak_list)
	# do peak detection on pre-trimmed data
	# trim by relative intensity
	apl = rel_threshold(peak_list, 2, copy_peaks=False)
	
	# trim by threshold
	new_peak_list = num_ions_threshold(apl, 3, 3000, copy_peaks=False)
	
	# ignore TMS ions and set mass range
	for peak in new_peak_list:
		peak.crop_mass(50, 400)
		peak.null_mass(73)
		peak.null_mass(147)
	
		# find area
		area = peak_sum_area(im_i, peak)
		peak.area = area
		area_dict = peak_top_ion_areas(im_i, peak)
		peak.ion_areas = area_dict
	
	return new_peak_list


@pytest.fixture(scope="session")
def peak(im_i):
	scan_i = im_i.get_index_at_time(31.17 * 60.0)
	ms = im_i.get_ms_at_index(scan_i)
	return Peak(12.34, ms)


@pytest.fixture(scope="session")
def ms(im_i):
	return deepcopy(im_i.get_ms_at_index(0))
	
	
@pytest.fixture(scope="function")
def expr(filtered_peak_list):
	# create an experiment
	return Experiment("ELEY_1_SUBTRACT", filtered_peak_list)


