"""proc.py
"""

from pyms.GCMS.IO.JCAMP import JCAMP_reader
from pyms.GCMS.Function import build_intensity_matrix

# read the raw data as a GCMS_data object
jcamp_file = "data/gc01_0812_066.jdx"
data = JCAMP_reader(jcamp_file)

# IntensityMatrix
# must build intensity matrix before accessing any intensity matrix methods.

# bin interval of 0.5, eg. for double charge ions
print("intensity matrix, bin interval = 0.5, boundary +/- 0.25")
im = build_intensity_matrix(data, 0.5, 0.25, 0.25)

print(" -> size of intensity matrix (#scans, #bins):", im.get_size())

print(" -> start mass:", im.min_mass)
print(" -> end mass:", im.max_mass)

index = im.get_index_of_mass(73.3)
print(" -> the index of the nearest mass to 73.3m/z is:", index)
print(" -> the nearest mass to 73.3m/z is:", im.get_mass_at_index(index))

