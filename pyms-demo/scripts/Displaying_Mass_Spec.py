"""
proc.py

Plot Mass Spectrum using matplotlib
"""

import sys
sys.path.append("../..")

import pathlib

import matplotlib.pyplot as plt

from pyms.GCMS.IO.JCAMP import JCAMP_reader
from pyms.IntensityMatrix import build_intensity_matrix_i
from pyms.Display import plot_mass_spec


data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
# Change this if the data files are stored in a different location

output_directory = pathlib.Path(".").resolve() / "output"

jcamp_file = data_directory / "gc01_0812_066.jdx"
data = JCAMP_reader(jcamp_file)
tic = data.tic
im = build_intensity_matrix_i(data)

# Extract the desired mass spectrum
ms = im.get_ms_at_index(1024)

fig, ax = plt.subplots(1, 1)

# Plot the spectrum
plot_mass_spec(ax, ms)

# Set the title
ax.set_title("Mass Spectrum at index 1024")

# Reduce the x-axis range to better visualise the data
ax.set_xlim(50, 350)

plt.show()
