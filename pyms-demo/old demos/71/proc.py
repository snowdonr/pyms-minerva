"""
proc.py

Plot detected peaks using matplotlib
"""

import sys
sys.path.append("../..")

import pathlib

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

from pyms.BillerBiemann import BillerBiemann, rel_threshold, num_ions_threshold
from pyms.Display import plot_ic, plot_peaks
from pyms.GCMS.IO.JCAMP import JCAMP_reader
from pyms.IntensityMatrix import build_intensity_matrix
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.TopHat import tophat


data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
# Change this if the data files are stored in a different location

output_directory = pathlib.Path(".").resolve() / "output"

jcamp_file = data_directory / "gc01_0812_066.jdx"
data = JCAMP_reader(jcamp_file)
data.trim("500s", "2000s")
tic = data.tic
im = build_intensity_matrix(data)

# Perform pre-filtering and peak detection.

n_scan, n_mz = im.size

for ii in range(n_mz):
	ic = im.get_ic_at_index(ii)
	ic_smooth = savitzky_golay(ic)
	ic_bc = tophat(ic_smooth, struct="1.5m")
	im.set_ic_at_index(ii, ic_bc)

# Detect Peaks
peak_list = BillerBiemann(im, points=9, scans=2)

print("Number of peaks found: ", len(peak_list))

# Filter the peak list, first by removing all intensities in a peak less than a
# given relative threshold, then by removing all peaks that have less than a
# given number of ions above a given value

pl = rel_threshold(peak_list, percent=2)
new_peak_list = num_ions_threshold(pl, n=3, cutoff=10000)

print("Number of filtered peaks: ", len(new_peak_list))

# Get Ion Chromatograms for 4 separate m/z channels.
ic191 = im.get_ic_at_mass(191)
ic73 = im.get_ic_at_mass(73)
ic57 = im.get_ic_at_mass(57)
ic55 = im.get_ic_at_mass(55)

# Create a subplot and plot the TIC and peaks
ms = im.get_ms_at_index(1024)

fig, ax = plt.subplots(1, 1)

# Plot the ICs
plot_ic(ax, tic, label="TIC")
plot_ic(ax, ic191, label="m/z 191")
plot_ic(ax, ic73, label="m/z 73")
plot_ic(ax, ic57, label="m/z 57")
plot_ic(ax, ic55, label="m/z 55")

# Plot the peaks
plot_peaks(ax, new_peak_list)

# Set the title
ax.set_title('TIC, ICs, and PyMS Detected Peaks')

# Add the legend
plt.legend()

plt.show()
