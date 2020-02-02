#################################################################
#       Creates Peak list for use with Display module
#################################################################

import pathlib

from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.IntensityMatrix import build_intensity_matrix_i
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.TopHat import tophat
from pyms.Display import Display
from pyms.BillerBiemann import BillerBiemann, rel_threshold, num_ions_threshold

data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
# Change this if the data files are stored in a different location

output_directory = pathlib.Path(".").resolve() / "output"

# read raw data
andi_file = data_directory / "MM-10.0_1_no_processing.cdf"
data = ANDI_reader(andi_file)

# Build Intensity Matrix
im = build_intensity_matrix_i(data)
n_scan, n_mz = im.size

# perform necessary pre filtering
for ii in range(n_mz):
    ic = im.get_ic_at_index(ii)
    ic_smooth = savitzky_golay(ic)
    ic_bc = tophat(ic_smooth, struct="1.5m")
    im.set_ic_at_index(ii, ic_bc)
    
# Detect Peaks
peak_list = BillerBiemann(im, points=3, scans=2)
print("Number of peaks found: ", len(peak_list))

######### Filter peaks###############
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

# TIC from raw data
tic = data.get_tic()
# save TIC to a file

# Get Ion Chromatograms for all m/z channels
n_mz = len(im.get_mass_list())
ic_list = []

for m in range(n_mz):
    ic_list.append(im.get_ic_at_index(m))

# Create a new display object, this time plot the ICs 
# and the TIC, as well as the peak list
display = Display()
display.plot_tic(tic, 'TIC')
for ic in ic_list:
    display.plot_ic(ic)
display.plot_peaks(new_peak_list, 'Peaks')
display.do_plotting('TIC, and PyMassSpec Detected Peaks')
display.show_chart()
