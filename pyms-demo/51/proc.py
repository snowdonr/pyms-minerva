"""proc.py
"""

from pyms.GCMS.Function import build_intensity_matrix_i
from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.Peak.Class import Peak

# read file and convert to intensity matrix
andi_file = "data/gc01_0812_066.cdf"
data = ANDI_reader(andi_file)
im = build_intensity_matrix_i(data)

# Get the scan of a known TIC peak (at RT 31.17 minutes)
# get the index of the scan nearest to 31.17 minutes (converted to seconds)
scan_i = im.get_index_at_time(31.17*60.0)
# get the MassSpectrum Object
ms = im.get_ms_at_index(scan_i)

# create a Peak object
peak = Peak(31.17, ms, minutes=True)

print(peak.get_UID())

# modify the range and null TMS ions
peak.crop_mass(60, 450)
peak.null_mass(73)
peak.null_mass(147)

# New UID after modification
print(peak.get_UID())

# Create another peak from an isomer of the first peak (at RT 31.44 minutes)
scan_i = im.get_index_at_time(31.44*60.0)
ms = im.get_ms_at_index(scan_i)

peak2 = Peak(31.44, ms, minutes=True)
# modify the range and null TMS ions
peak2.crop_mass(60, 450)
peak2.null_mass(73)
peak2.null_mass(147)

# Compare isomers UID
print(peak2.get_UID())
