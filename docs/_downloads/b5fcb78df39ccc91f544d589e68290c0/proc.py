"""proc.py
"""

from pyms.GCMS.IO.ANDI.Function import ANDI_reader
from pyms.GCMS.Function import build_intensity_matrix_i

# read the raw data
andi_file = "data/gc01_0812_066.cdf"

data = ANDI_reader(andi_file)

# original mass range
print(data.get_min_mass(), data.get_max_mass())

# Build an integer mass intensity matrix
im = build_intensity_matrix_i(data)

# Crop the mass range to 60 to 400 m/z
im.crop_mass(60, 400)

# new mass range
print(im.get_min_mass(), im.get_max_mass())

# Set intensities of masses related to solvents to null (zero)
im.null_mass(73)
im.null_mass(147)

# test if all intensities zeroed
print(sum(im.get_ic_at_mass(73).get_intensity_array()))
print(sum(im.get_ic_at_mass(147).get_intensity_array()))
