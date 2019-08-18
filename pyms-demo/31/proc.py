"""proc.py
"""

from pyms.GCMS.IO.JCAMP.Function import JCAMP_reader
from pyms.GCMS.Function import build_intensity_matrix
#from pyms.Utils.IO import save_data

# read the raw data as a GCMS_data object
jcamp_file = "data/gc01_0812_066.jdx"
data = JCAMP_reader(jcamp_file)

# IntensityMatrix
# must build intensity matrix before accessing any intensity matrix methods.

# default, float masses with interval (bin interval) of one from min mass
print("default intensity matrix, bin interval = 1, boundary +/- 0.5")
im = build_intensity_matrix(data)

#
# MassSpectrum
#

ms = im.get_ms_at_index(0)

# attributes and properties
print(len(ms))
print(len(ms.mass_list))
print(len(ms.mass_spec))

#
# IonChromatogram
#

# TIC from raw data
tic = data.get_tic()
# save TIC to a file
tic.write("output/tic.dat",minutes=True)

# get the first ion chromatogram of the IntensityMatrix
ic = im.get_ic_at_index(0)
ic.write("output/ic_index_0.dat",minutes=True)
# get the ion chromatogram for m/z = 73
ic = im.get_ic_at_mass(73)
ic.write("output/ic_mass_73.dat",minutes=True)

# some tests on ion chromatogram objects
print("'tic' is a TIC:", tic.is_tic())
print("'ic' is a TIC:", ic.is_tic())

