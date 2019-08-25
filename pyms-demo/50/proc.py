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

# Get the retention time (in seconds)
print(peak.rt)

# Get the peaks unique ID
# Consists of the two most abundant ions and their ratio,
# and the retention time (in the format set by minutes=True or False)
print(peak.UID)

# Create another peak from an isomer of the first peak (at RT 31.44 minutes)
scan_i = im.get_index_at_time(31.44*60.0)
ms = im.get_ms_at_index(scan_i)
peak2 = Peak(31.44, ms, minutes=True)
print(peak2.UID)
