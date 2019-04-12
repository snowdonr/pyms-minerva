"""proc.py
"""

from pyms.GCMS.IO.JCAMP.Function import JCAMP_reader
from pyms.GCMS.Function import build_intensity_matrix

# read the raw data as a GCMS_data object
jcamp_file = "data/gc01_0812_066.jdx"
data = JCAMP_reader(jcamp_file)

# IntensityMatrix
# must build intensity matrix before accessing any intensity matrix methods.

# default, float masses with interval (bin interval) of one from min mass
print("default intensity matrix, bin interval = 1, boundary +/- 0.5")
im = build_intensity_matrix(data)

print("size of intensity matrix (#scans, #bins):", im.get_size())

print("start mass:", im.get_min_mass())
print("end mass:", im.get_max_mass())

index = im.get_index_of_mass(73.3)
print("the index of the nearest mass to 73.3m/z is:", index)
print("the nearest mass to 73.3m/z is:", im.get_mass_at_index(index))

# get the list of masses (bin centers), and print the first ten
masses = im.get_mass_list()
print(masses[:10])

