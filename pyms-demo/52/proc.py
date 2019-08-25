"""proc.py
"""

from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.GCMS.Function import build_intensity_matrix_i
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.Baseline.TopHat import tophat

from pyms.Deconvolution.BillerBiemann.Function import BillerBiemann

# read the raw data as a GCMS_data object
andi_file = "data/gc01_0812_066.cdf"
data = ANDI_reader(andi_file)

im = build_intensity_matrix_i(data)

n_scan, n_mz = im.size

print("Intensity matrix size (scans, masses):", (n_scan, n_mz))

# noise filter and baseline correct
for ii in range(n_mz):
    ic = im.get_ic_at_index(ii)
    ic_smooth = savitzky_golay(ic)
    ic_bc = tophat(ic_smooth, struct="1.5m")
    im.set_ic_at_index(ii, ic_bc)

# Use Biller and Biemann technique to find apexing ions at a scan
# default is maxima over three scans and not to combine with any neighbouring
# scan.
peak_list = BillerBiemann(im)

print("Number of peaks found: ", len(peak_list))


# Find apex oven 9 points and combine with neighbouring peak if two scans apex
# next to each other.
peak_list = BillerBiemann(im, points=9, scans=2)

print("Number of peaks found: ", len(peak_list))
