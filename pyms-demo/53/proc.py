"""proc.py
"""

from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.GCMS.Function import build_intensity_matrix_i
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.Baseline.TopHat import tophat

from pyms.Deconvolution.BillerBiemann.Function import BillerBiemann, \
    rel_threshold, num_ions_threshold

# read the raw data as a GCMS_data object
andi_file = "data/gc01_0812_066.cdf"
data = ANDI_reader(andi_file)

im = build_intensity_matrix_i(data)

n_scan, n_mz = im.get_size()

print("Intensity matrix size (scans, masses):", (n_scan, n_mz))

# noise filter and baseline correct
for ii in range(n_mz):
    ic = im.get_ic_at_index(ii)
    ic_smooth = savitzky_golay(ic)
    ic_bc = tophat(ic_smooth, struct="1.5m")
    im.set_ic_at_index(ii, ic_bc)

# Use Biller and Biemann technique to find apexing ions at a scan.
# Find apex oven 9 points and combine with neighbouring peak if two scans apex
# next to each other.
peak_list = BillerBiemann(im, points=9, scans=2)

print("Number of peaks found: ", len(peak_list))

# Filter the peak list,
# first by removing all intensities in a peak less than a given relative
# threshold,
# then by removing all peaks that have less than a given number of ions above
# a given value

# Parameters
# percentage ratio of ion intensity to max ion intensity
r = 2

# minimum number of ions, n
n = 3
# greater than or equal to threshold, t
t = 10000

# trim by relative intensity
pl = rel_threshold(peak_list, r)

# trim by threshold
new_peak_list = num_ions_threshold(pl, n, t)

print("Number of filtered peaks: ", len(new_peak_list))
